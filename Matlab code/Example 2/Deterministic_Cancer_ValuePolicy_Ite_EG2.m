% This function computes the value function of the deterministic
% optimization problem of adaptive cancer therapy of our 2nd example
% Carrère 2017 https://doi.org/10.1016/j.jtbi.2016.11.009
% in the following form:
% u = max_d{tau*(dmax + delta) + u_interp} over d in {0,dmax})
% We recover the optimal bang-bang policy in feedback form via a
% semi-Lagrangian discretization.
% This function applies Gauss-Seidel iterations
% and applies a "value-policy" iteration scheme described in the Supplementary Materials
% of the paper "Threshold-awareness in adaptive cancer therapy"
% (https://www.biorxiv.org/content/10.1101/2022.06.17.496649v2.supplementary-material)
%
% N (input): number of points (minus 1) along each side of the unit qp-square
%
% U(output): a (N+1)x(N+1) matrix of the value function
% Dmat (output): a (N+1)x(N+1) matrix of the optimal policy in feedback form
%
% Author: MingYi Wang, Cornell University
% Last Modified: 12/2023
%
function  [U,Dmat] = Deterministic_Cancer_ValuePolicy_Ite_EG2(N)
%% parameters
dx = 1/N; %spatial discretization along one side of qp-square
dy = dx; %uniform spatial grid
h = sqrt(dx^2+dy^2);
largeNum=10^4; %use a large number to substitue the "infinity" BC in practice
tol= 10^-6; %tolerance of convergence
xx = 0:dx:1;
%parameters from Carrere 2017 https://doi.org/10.1016/j.jtbi.2016.11.009
time_factor = 1;
gs = 0.031*time_factor;
gr = 0.026*time_factor;
m = 30;
K = 0.48; %of order 1e7
beta = 6.25*time_factor; %of order 1e-7
beta_hat = K*beta;
alpha = 0.06*time_factor;

dmax = 3;
s = 0.05;

r_b = 0.01;
f_b = 1- r_b;
ind_rec=length(0:dx:r_b);%index for the recovery barrier
ind_death=find(xx==f_b,1);%index for the death barrier

%drift functions
fq = @(q,p,d) (1-p).*q.*(1-q).*(gs-gr) + beta_hat.*p.*q.^2.*(1-q) - alpha.*q.*(1-q).*d;
fp = @(q,p,d) p.*(1-p).*(gs.*q + gr.*(1-q)) - alpha.*q.*p.*d - beta_hat.*p.^2.*q.*(1-q);
%preset factor, if the infinity-norm between two value iterations falls
%below the previously stored error record, we enter the "policy-update" step
rho = 0.7;

%time step to compute the foot of characteristics
tau = 0.15*sqrt(h);
%% Initialization
U = zeros(N+1,N+1);
%Boundary conditions
U(1:ind_rec-1,:) = 0;
U(ind_death+1:end,:) = largeNum;
U(ind_rec:end,1) = largeNum;
%initialize the interior domain to some arbitrary poistive values
U(ind_rec:ind_death,2:N+1) = largeNum + 1;

Dmat = zeros(N+1,N+1); %initialize the policy matrix to be the zero matrix

%% precompute coefficients in interior domain

tic
[XQ1,XQ2,XP1,XP2,KQ1,KQ2,KP1,KP2]...
    = precompute_coeff_2D(fp,fq,dmax,N,r_b,f_b,ind_rec,ind_death,tau);
t1 = toc
%% fast sweeping method
%------------------------Initialization------------------------------------
k=0; %value iteration counter
kk=0; %"policy-update" counter
err=100; %initialize an error of convergence
Inf_err_list = [];
PE_flag = false;
%----------------------------Main Loop---------------------------------
tic
iOrder=(ind_rec):1:(ind_death);
jOrder=N+1:-1:2;
while err>tol
    uold=U;
    for i = iOrder
        for j = jOrder
            %Since we have the bang-bang property, we only have d=0 or
            %d=dmax, exclusively.
            %First, we consider the case when d = 0 (no therapy)
            x1q=XQ1(i,j);
            x1p=XP1(i,j);
            kq=KQ1(i,j);
            kp=KP1(i,j);

            if x1p < r_b
                uij1 = 0;
            elseif x1p > f_b
                uij1 = largeNum;
            else

                ql= xx(kq-1);qr=xx(kq);pb=xx(kp-1);pt=xx(kp);
                u11=U(kp-1,kq-1);u12=U(kp-1,kq);
                u21=U(kp,kq-1);u22=U(kp,kq);
                uij1 = BilinearInterp(u11,u12,u21,u22,ql,qr,pt,pb,x1q,x1p,dx,dy);


            end
            u1 = tau*s+uij1; %running cost + interpolated value


            %Next, we consider the case when d = dmax (MTD rate)
            x2q=XQ2(i,j);
            x2p=XP2(i,j);
            kq=KQ2(i,j);
            kp=KP2(i,j);

            if x2p < r_b
                uij2 = 0;
            elseif x2p > f_b
                uij2 = largeNum;
            else


                ql=xx(kq-1);qr=xx(kq);pb=xx(kp-1);pt=xx(kp);
                u11=U(kp-1,kq-1);u12=U(kp-1,kq);
                u21=U(kp,kq-1);u22=U(kp,kq);
                uij2 = BilinearInterp(u11,u12,u21,u22,ql,qr,pt,pb,x2q,x2p,dx,dy);


            end
            u2 = tau*(dmax+s)+uij2; %running cost + interpolated value



            unew=min(u1,u2); %take the smaller value
            if u1 <= u2
                d_new = 0;
            else
                d_new = 1;
            end
            if unew < U(i,j)
                U(i,j)=unew;
                Dmat(i,j) = d_new;
            end
        end
    end

    if k == 10
        %store the initial error after 4 sweeps
        errlist(kk+1) = max(max(abs(uold-U)));
        err = errlist(kk+1);
    end
    %Printing
    k
    Inf_Norm_per_ite = max(max(abs(uold-U)))
    Inf_err_list = [Inf_err_list,Inf_Norm_per_ite];

    if ((k > 10) && (Inf_Norm_per_ite < 200))
        diff = Inf_Norm_per_ite;
        if (diff < rho*err)
            PE_flag = true;

        else
            PE_flag = false;
        end

        if PE_flag %criterion to enter the "policy-update" step
            kk = kk+1
            err = min(err,diff);
            errlist(kk+1) = err;
            

            %--------------------"policy-evaluation" step--------------------------
            u_vec = U(:); %convert the value function to be a long vector
            d_vec = Dmat(:); %convert the policy matrix to be a long vector
            b_vec = zeros(1,(N+1)*(N+1)); %initialize b vector

            [column,row]=meshgrid(1:N+1,1:N+1);
            xloc_mat = zeros(N+1,N+1);
            yloc_mat = zeros(N+1,N+1);
            kx_mat= zeros(N+1,N+1);
            ky_mat = zeros(N+1,N+1);
            tau_mat = zeros(N+1,N+1);

            %first, take care of the ones on the diagonal
            I = 1:(N+1)*(N+1);
            J = 1:(N+1)*(N+1);
            V = ones(1, (N+1)*(N+1));

            %finding nodes where d=0 or d=dmax
            test_mtd = (Dmat==1);
            test_0 = (Dmat==0);

            %find locations and coefficients for d=dmax case
            xloc_mat(test_mtd) = XQ2(test_mtd);
            yloc_mat(test_mtd) = XP2(test_mtd);
            kx_mat(test_mtd)= KQ2(test_mtd);
            ky_mat(test_mtd) = KP2(test_mtd);


            %find locations and coefficients for d=0 case
            xloc_mat(test_0) = XQ1(test_0);
            yloc_mat(test_0) = XP1(test_0);
            kx_mat(test_0)= KQ1(test_0);
            ky_mat(test_0) = KP1(test_0);

            %Taking care of Boundary conditions
            test_u0 = (U == 0);
            node_index = (reshape(column(test_u0), 1,[])-1)*(N+1) + reshape(row(test_u0), 1,[]);
            b_vec(node_index) = u_vec(node_index);

            test_largenum = (U == largeNum);
            node_index = (reshape(column(test_largenum), 1,[])-1)*(N+1) + reshape(row(test_largenum), 1,[]);
            b_vec(node_index) = u_vec(node_index);

            test_ky0 = (ky_mat == 0);
            node_index = (reshape(column(test_ky0), 1,[])-1)*(N+1) + reshape(row(test_ky0), 1,[]);
            b_vec(node_index) = u_vec(node_index);

            t_int = ~(test_u0 | test_largenum | test_ky0);

            %Taking care of the rest
            t = (t_int);
            node_index = (reshape(column(t), 1,[])-1)*(N+1) + reshape(row(t), 1,[]);
            xloc = xloc_mat(t)';
            yloc = yloc_mat(t)';
            kx = kx_mat(t);
            ky = ky_mat(t);



            xl=xx(kx-1);
            xr=xx(kx);
            yb=xx(ky-1);
            yt=xx(ky);
            col_indx1 = (kx-2)*(N+1)+(ky-1);
            col_indx2 = col_indx1+(N+1);
            col_indx3 = col_indx1+1;
            col_indx4 = col_indx3+(N+1);
            %weights of bilinear interpolation in 2D
            weight1 = -((yt-yloc)/dy).*((xr-xloc)/dx);
            weight2 = -((yt-yloc)/dy).*((xloc-xl)/dx);
            weight3 = -((yloc-yb)/dy).*((xr-xloc)/dx);
            weight4 = -((yloc-yb)/dy).*((xloc-xl)/dx);
            I = [I node_index node_index node_index node_index];
            J = [J col_indx1' col_indx2' col_indx3' col_indx4'];
            V = [V weight1 weight2 weight3 weight4];
            b_vec(node_index) = tau.*(dmax.*d_vec(node_index)+s);

            A = sparse(I,J,V,(N+1)*(N+1),(N+1)*(N+1));

            u_true = A\b_vec';
            U = reshape(u_true,[N+1 N+1]);
        end
    end
    k=k+1;
    err
end
t=toc

%% Visualization
[X,Y] = meshgrid(xx,xx);
UU = U;
UU(U>= largeNum) = NaN;
figure
contourf(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),UU(1:2:end,1:2:end),15,'EdgeColor','none');
hold on
yline(r_b,'c:','LineWidth',1.5);
yline(f_b,'c:','LineWidth',1.5);
hold off
xlabel('Fraction of the sensitive (q)','FontSize',14);
ylabel('Normalized total size (p)','FontSize',14);
colorbar();

figure
contourf(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),Dmat(1:2:end,1:2:end),[0,1]);
hold on
yline(r_b,'c:','LineWidth',1.5);
yline(f_b,'c:','LineWidth',1.5);
hold off
xlabel('Fraction of the sensitive (q)','FontSize',14);
ylabel('Normalized total size (p)','FontSize',14);

end
%% subfunction
%This function precomputes coefficients in interior domain
%
% fq,fq (input): anonymous functions of the deterministic cancer dynamics
% dmax (input): MTD rate
% N (input): number of sample points along each side of the unit square
% ind_rec (input): row index of the recovery barrier
% ind_death (input): row index of the death barrier
% tau (input): time step to compute the foot of characteristics
%
% XQ1,XQ2,XP1,XP2 (output): matrices of foot of characteristics in q- and p-
%                          direction respectively.
% KQ1,KQ2,KP1,KP2 (output): coefficients of foot of characteristics in
%                           q- and p-direction respectively.
% Note: For the output above: "1" for d=0 while 2 for d=dmax.
%
function [XQ1,XQ2,XP1,XP2,KQ1,KQ2,KP1,KP2]...
    =precompute_coeff_2D(fp,fq,dmax,N,r_b,f_b,ind_rec,ind_death,tau)
%initialization of storage
XQ1=zeros(N+1,N+1);
XQ2=zeros(N+1,N+1);
XP1=zeros(N+1,N+1);
XP2=zeros(N+1,N+1);
KP1=zeros(N+1,N+1);
KP2=zeros(N+1,N+1);
KQ1=zeros(N+1,N+1);
KQ2=zeros(N+1,N+1);

%create a uniform grid on [0,1]x[0,1]
q = linspace(0,1,N+1);
p = q;
for i=(ind_rec):(ind_death)
    for j=2:N+1
        %Since we have bang-bang property, we only consider two
        %controls: d=0 or d=dmax=3 exclusively

        %We first compute the foot of characteristics for d=0
        fx1 = fq(q(j),p(i),0);
        fy1 = fp(q(j),p(i),0);
        q_ij1 = q(j)+tau*fx1;
        p_ij1 = p(i)+tau*fy1;
        if p_ij1 >= r_b && p_ij1 <= f_b

            %compute the coefficients
            kq = find(q_ij1<=q',1);
            kp = find(p_ij1<=p',1);
            %store them in the corresponding matrices
            KQ1(i,j) = kq;
            KP1(i,j) = kp;
        end
        XQ1(i,j) = q_ij1;
        XP1(i,j) = p_ij1;


        %Repeat the process for d = damx
        fx2 = fq(q(j),p(i),dmax);
        fy2 = fp(q(j),p(i),dmax);
        q_ij2 = q(j)+tau*fx2;
        p_ij2 = p(i)+tau*fy2;

        if p_ij2 >= r_b && p_ij2 <= f_b
            kq = find(q_ij2<=q',1);
            kp = find(p_ij2<=p',1);

            KQ2(i,j)=kq;
            KP2(i,j)=kp;

        end
        XQ2(i,j)=q_ij2;
        XP2(i,j)=p_ij2;
    end
end
end




