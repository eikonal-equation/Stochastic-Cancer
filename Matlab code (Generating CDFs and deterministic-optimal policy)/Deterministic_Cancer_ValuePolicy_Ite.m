% This function computes the value function of the deterministic
% optimization problem of adaptive cancer therapy from Gluzman et al.
% https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
% in the following form:
% u = max_d{tau*(dmax + delta) + u_interp} over d in {0,dmax}
% We hence recover the optimal bang-bang policy in feedback form via a
% semi-Lagrangian discretization.
% This function appliesGauss-Seidel iterations with "Fast sweeping" ordering
% and applies a "value-policy" iteration scheme described in our paper
% (paper name & link)
%
% N (input): number of points (minus 1) along each side of the unit qp-square
%
% U(output): a (N+1)x(N+1) matrix of the value function
% Dmat (output): a (N+1)x(N+1) matrix of the optimal policy in feedback form
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/18/2022
%
function  [U,Dmat] = Deterministic_Cancer_ValuePolicy_Ite(N)
%% parameters
dx = 1/N; %spatial discretization along one side of qp-square
dy = dx; %uniform spatial grid
h = sqrt(dx^2+dy^2);
largeNum=10^5; %use a large number to substitue the "infinity" BC in practice
tol= 10^-6; %tolerance of convergence
xx = 0:dx:1;
% parameters from Gluzman et al. https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
ba = 2.5;bv = 2; c = 1; n = 4;
s=0.05;%treatment cost
dmax=3;%MTD rate

ind_rec=length(0:dx:0.01);%index for the recovery barrier
ind_death=find(xx==0.99,1);%index for the death barrier

%drift functions in q(horizontal)- and p(vertical)-directions.
fp=@(q,p,d) p*(1-p)*(ba/(n+1)-q*(bv-c)-d);
fq=@(q,p) q*(1-q)*(bv/(n+1)*(1+p+p^2+p^3+p^4)-c);

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
%initialize the interior domain to some arbitrary poistive values
U(ind_rec:ind_death,:) = 20;

Dmat = zeros(N+1,N+1); %initialize the policy matrix to be the zero matrix

%% precompute coefficients in interior domain

tic
[XQ1,XQ2,XP1,XP2,KQ1,KQ2,KP1,KP2]...
    = precompute_coeff_2D(fp,fq,dmax,N,ind_rec,ind_death,tau);
t1 = toc
%% fast sweeping method
%------------------------Initialization------------------------------------
k=0; %value iteration counter
kk=0; %"policy-update" counter
err=100; %initialize an error of convergence

%----------------------------Main Loop---------------------------------
tic
while err>tol
    uold=U;
    switch mod(k,4) %sweep in 4 directions
        case 0
            iOrder=(ind_rec):1:(ind_death);
            jOrder=1:1:N+1;
        case 1
            iOrder=(ind_rec):1:(ind_death);
            jOrder=N+1:-1:1;
        case 2
            iOrder=(ind_death):-1:(ind_rec);
            jOrder=N+1:-1:1;
        case 3
            iOrder=(ind_death):-1:(ind_rec);
            jOrder=1:1:N+1;
    end
    
    for i = iOrder
        for j = jOrder
            %Since we have the bang-bang property, we only have d=0 or
            %d=dmax, exclusively.
            %First, we consider the case when d = 0 (no therapy)
            x1q=XQ1(i,j);
            x1p=XP1(i,j);
            kq=KQ1(i,j);
            kp=KP1(i,j);
            
            if j==1 % on q = 0
                uij1 = LinearInterp(U(:,1),kp,x1p,dy,0);
            elseif j==N+1 %on q = 1
                uij1= LinearInterp(U(:,N+1),kp,x1p,dy,0);
                
            else %in the interior
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
            
            
            if j==1 %on q = 0
                uij2 = LinearInterp(U(:,1),kp,x2p,dy,0);
            elseif j==N+1 %on q = 1
                uij2= LinearInterp(U(:,N+1),kp,x2p,dy,0);
                
            else %in the interior
                ql=xx(kq-1);qr=xx(kq);pb=xx(kp-1);pt=xx(kp);
                u11=U(kp-1,kq-1);u12=U(kp-1,kq);
                u21=U(kp,kq-1);u22=U(kp,kq);
                uij2 = BilinearInterp(u11,u12,u21,u22,ql,qr,pt,pb,x2q,x2p,dx,dy);
                
            end
            u2 = tau*(dmax+s)+uij2; %running cost + interpolated value
            
            unew=min(u1,u2); %take the smaller value
            if u1<=u2
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
    if k==3
        %store the initial error after 4 sweeps
        errlist(kk+1) = max(max(abs(uold-U)));
        err = errlist(kk+1);
    end
    
    k
    Inf_Norm_per_ite = max(max(abs(uold-U)))
    if (k > 3) && (mod(k,4) == 3)
        diff = max(max(abs(uold-U)));
        if diff < rho*err %criterion to enter the "policy-update" step
            kk = kk+1
            errlist(kk+1)=diff;
            err = diff;
            
            %------------------"policy-update"----------------------------
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
            test_mtd = (Dmat==dmax);
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
            test_u0 = (U==0);
            node_index = (reshape(column(test_u0), 1,[])-1)*(N+1) + reshape(row(test_u0), 1,[]);
            b_vec(node_index) = u_vec(node_index);
            
            test_largenum = (U>=largeNum);
            node_index = (reshape(column(test_largenum), 1,[])-1)*(N+1) + reshape(row(test_largenum), 1,[]);
            b_vec(node_index) = u_vec(node_index);
            
            test_ky0 = (ky_mat == 0);
            node_index = (reshape(column(test_ky0), 1,[])-1)*(N+1) + reshape(row(test_ky0), 1,[]);
            b_vec(node_index) = u_vec(node_index);
            
            t_int = ~(test_u0 | test_largenum | test_ky0);
            
            %Taking care the case on q = 0
            t = (t_int & (xloc_mat == 0));
            node_index = (reshape(column(t), 1,[])-1)*(N+1) + reshape(row(t), 1,[]);
            yloc = yloc_mat(t)';
            ky = ky_mat(t);
            yb = xx(ky-1);
            yt = xx(ky);
            
            
            
            col_indx1 = ky-1;
            col_indx2 = ky;
            weight1 = -(yt-yloc)/dy;%weights of linear interpolation in 1D
            weight2 = -(yloc-yb)/dy;
            I = [I node_index node_index];
            J = [J col_indx1' col_indx2'];
            V = [V weight1 weight2];
            b_vec(node_index) = tau.*(d_vec(node_index)+s);
            
            
            %Taking care the case on q = 1
            t = (t_int & (xloc_mat == 1));
            node_index = (reshape(column(t), 1,[])-1)*(N+1) + reshape(row(t), 1,[]);
            yloc = yloc_mat(t)';
            ky = ky_mat(t);
            yb = xx(ky-1);
            yt = xx(ky);
            
            
            col_indx1 = N*(N+1)+(ky-1);
            col_indx2 = N*(N+1)+ky;
            weight1 = -(yt-yloc)/dy;%weights of linear interpolation in 1D
            weight2 = -(yloc-yb)/dy;
            I = [I node_index node_index];
            J = [J col_indx1' col_indx2'];
            V = [V weight1 weight2];
            b_vec(node_index) = tau.*(d_vec(node_index)+s);
            
            
            %Taking care of the interior
            t = (t_int &  ~((xloc_mat==0) | (xloc_mat==1)));
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
            b_vec(node_index) = tau.*(d_vec(node_index)+s);
            
            A = sparse(I,J,V,(N+1)*(N+1),(N+1)*(N+1));
            
            u_true = A\b_vec';
            U = reshape(u_true,[N+1 N+1]);
        end
    end
    k=k+1;
    err
end
t=toc
%% Visualization of the optimal policy
if N >= 800 %if N >= 800, we suggest downsampling it to 801x801
    NN = 800;
    q=linspace(0,1,NN+1);
    p=linspace(0,1,NN+1);
    
    %covert the qp-square into GLY-DEF-VOP triangle
    x1=[];x2=[];x3=[];
    for ii=1:1:NN+1
        x1ij=p(ii)*ones(1,NN+1);
        x2ij=(1-p(ii))*(1-q);
        x3ij=1-x1ij-x2ij;
        x1=[x1 x1ij];
        x2=[x2 x2ij];
        x3=[x3 x3ij];
    end
    [x, y] = terncoords(x3, x2, x1);
    
    
    Xtri=zeros(NN+1,NN+1); Ytri=zeros(NN+1,NN+1);
    for jj=1:NN+1
        Xtri(jj,1:NN+1)=x(NN*(jj-1)+jj:NN*(jj-1)+jj+NN);
        Ytri(jj,1:NN+1)=y(NN*(jj-1)+jj:NN*(jj-1)+jj+NN);
    end
else
    q=linspace(0,1,N+1);
    p=linspace(0,1,N+1);
    %covert the qp-square into GLY-DEF-VOP triangle
    x1=[];x2=[];x3=[];
    for ii=1:1:N+1
        x1ij=p(ii)*ones(1,N+1);
        x2ij=(1-p(ii))*(1-q);
        x3ij=1-x1ij-x2ij;
        x1=[x1 x1ij];
        x2=[x2 x2ij];
        x3=[x3 x3ij];
    end
    [x, y] = terncoords(x3, x2, x1);
    
    
    Xtri=zeros(N+1,N+1); Ytri=zeros(N+1,N+1);
    for jj=1:N+1
        Xtri(jj,1:N+1)=x(N*(jj-1)+jj:N*(jj-1)+jj+N);
        Ytri(jj,1:N+1)=y(N*(jj-1)+jj:N*(jj-1)+jj+N);
    end
end


%--------------------------show plot--------------------------------
N_plot = size(Xtri,1)-1;
figure
[AA,BB]=contourf(Xtri,Ytri,Dmat(1:N/N_plot:end,1:N/N_plot:end),[0 1]);
set(BB,'LineColor','none');
hold on 
plot(linspace(0.4950,0.5050,3),[0.8574,0.8574,0.8574],'c:','linewidth',1.5);
plot(linspace(0.005,0.99,100),0.0087*ones(1,100),'c:','linewidth',1.5);
hold off
axis equal
axis([0 1 0 1])
text(-0.09,0,'VOP','FontSize',11);
text(0.53,0.85,'GLY','FontSize',11);
text(1.01,0,'DEF','FontSize',11);
titlename = sprintf('The Optimal DrugsOn-Off Region for N=%g',N);
title(titlename);
set(gca,'visible','off');
set(findall(gca,'type','text'),'visible','on');


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
            =precompute_coeff_2D(fp,fq,dmax,N,ind_rec,ind_death,tau)
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
            for j=1:N+1
                %Since we have bang-bang property, we only consider two
                %controls: d=0 or d=dmax=3 exclusively
                
                %We first compute the foot of characteristics for d=0
                fx1 = fq(q(j),p(i));
                fy1 = fp(q(j),p(i),0);
                q_ij1 = q(j)+tau*fx1;
                p_ij1 = p(i)+tau*fy1;
                
                %compute the coefficients
                kq = find(q_ij1<=q',1);
                kp = find(p_ij1<=p',1);
                %store them in the corresponding matrices
                KQ1(i,j) = kq;
                KP1(i,j) = kp;
                XQ1(i,j) = q_ij1;
                XP1(i,j) = p_ij1;
                
                
                %Repeat the process for d = damx
                fx2 = fx1; %since drugs don't affect the dynamics in q-direction
                fy2 = fp(q(j),p(i),dmax);
                q_ij2 = q(j)+tau*fx2;
                p_ij2 = p(i)+tau*fy2;
                kq = find(q_ij2<=q',1);
                kp = find(p_ij2<=p',1);
                
                KQ2(i,j)=kq;
                KP2(i,j)=kp;
                XQ2(i,j)=q_ij2;
                XP2(i,j)=p_ij2;
            end
        end
    end
end



