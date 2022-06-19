function Xcost = CDF_stationary_policy(Dmat_det,xloc,yloc,choice,sample_size)
%This function computes the CDF y = Pr(J <= s) that measures the
%probability of success where the cumulative cost J is within any positive
%threshold value s using the stationary policy computed from Gluzman et al.
%https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
%
%Dmat_det (input): the policy matrix of stationary policy computed from
%                  Gluzman et al.
%xloc (input): x (or q)-coordinate of the starting tumor configuration
%yloc (input): y (or p)-coordinate of the starting tumor configuration
%choice (input): choice of the determination strategy, only 3 options:
%                (i) 'conservative'; (ii) 'aggressive'; (iii) 'majority'.
%sample_size(input): sample size of Monte Carlo simulations
%
%Xcost (output): Array of the accumulative (random) cost of each sample path
%
%(We assume the initial policy at (xloc,yloc) is found before executing
%this function)
%
%% parameters
N = length(Dmat_det)-1; %number of points along one side of spatial grid
dx = 1/N;
% parameters from Gluzman et al.
% https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
ba=2.5;
bv = 2;
c = 1;
n = 4;

s1=0.15;s2=0.15;s3=0.15;%small diffusion constant
% s1=0.5;s2=0.5;s3=0.5;%large diffusion constant

s=0.05;%treatment cost
dmax=3;%MTD

xx = linspace(0,1,N+1);
ind_rec=length(0:dx:0.01);%index for recovery barrier
ind_death=find(xx==0.99,1);%index for death barrier

%drift functions (deterministic portion) of the cancer dynamics
fp = @(q,p,d) p*(1-p)*(ba/(n+1)-q*(bv-c)-d)-p*(s1^2*p-(s1^2*p^2+s2^2*(1-p)^2*(1-q)^2 ...
    + s3^2*(1-p)^2*q^2));
fq = @(q,p) q*(1-q)*(bv/(n+1)*(1+p+p^2+p^3+p^4)-c)+q*(1-q)*( (1-q)*s2^2-q*s3^2  );


%% Monte Carlo Simulations
%set the time step adaptively such that for each step, we spend 1
dt0 = 0.005/(s)/20; %1/20-th of 0.005 cost when not using drugs
dtmax = 0.005/(dmax+s)/10;%1/10-th of 0.005 cost

count_death=0; %initialize the number of deaths to 0

% tic
parfor ii = 1:sample_size
    cost=0; %accumulative cost
    ylist=[yloc];
    xlist=[xloc];
    w1=[0];w2=[0];w3=[0];%the 3D Brownian motion starts at (0,0,0)
    
    %initial policy (assuming known)
    policy = 0;
    %     policy = dmax
    
    dmax_set=[];
    d0_set=[1];
    policy_list=[policy];
    
    
    %     dt=dtmax;
    dt=dt0;
    
    j=1;
    while ylist(j) > 0.01
        j=j+1;
        %Sample Brownian motion as W_n - W_{n-1} = sqrt{dt}*N(0,1)
        w1(j)=w1(j-1)+sqrt(dt)*normrnd(0,1);
        w2(j)=w2(j-1)+sqrt(dt)*normrnd(0,1);
        w3(j)=w3(j-1)+sqrt(dt)*normrnd(0,1);
        %coordinates of the position from last step
        yval=ylist(j-1);
        xval=xlist(j-1);
        %Run EM step
        ylist(j) = yval+fp(xval,yval,policy)*dt+ yval*(1-yval)*s1*(w1(j)-w1(j-1))+...
            yval*(1-yval)*(1-xval)*s2*(w2(j)-w2(j-1))+...
            yval*(1-yval)*xval*s3*(w3(j)-w3(j-1));
        xlist(j) = xval+fq(xval,yval)*dt+xval*(1-xval)*(s3*(w3(j)-w3(j-1))-s2*(w2(j)-w2(j-1)));
        
        cost=cost+(policy+s)*dt; %accumulating cost
        
        %if we cross the failure barrier, stop
        if ylist(j) > 0.99
            count_death = count_death+1;
            %set the cost to be a large value (the theoretical value should be infinite)
            cost = 10^6;
            break
        end
        
        if ylist(j) >= 0.01
            %find the indices of the current position on the grid
            kq=find(xlist(j)<=xx',1);
            kp=find(ylist(j)<=xx',1);
            
            %----------Determination of the policy at next step----------------
            if (kp == (ind_death)) %if we are just below the death barrier
                d3 = Dmat_det(kp,kq-1);d4 = Dmat_det(kp,kq);
                switch choice
                    case 'conservative'
                        %super-conservative determination strategy
                        if ((d3+ d4) == 2)
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt = dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'aggressive'
                        %super-aggressive determination strategy
                        if ((d3 + d4) > 0)
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt = dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'majority'
                        %vote by majority determination strategy
                        if ((d3 + d4) > 0)
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt = dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                end
            else
                %if in the interior, find the values of feedback policies
                %surrounding the current position
                d1=Dmat_det(kp-1,kq-1);d2=Dmat_det(kp-1,kq);
                d3=Dmat_det(kp,kq-1);d4=Dmat_det(kp,kq);
                d_square=[d1,d2,d3,d4];
                
                switch choice
                    case 'conservative'
                        %super-conservative determination strategy
                        if sum(d_square)== 4
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt = dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'aggressive'
                        %super-aggressive determination strategy
                        if sum(d_square) > 0
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt=dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt=dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'majority'
                        %vote by majority determination strategy
                        if sum(d_square)>=3
                            policy = dmax;
                            
                            dmax_set=[dmax_set,j];
                            dt=dtmax;
                        else
                            policy = 0;
                            
                            d0_set=[d0_set,j];
                            dt=dt0;
                        end
                        policy_list=[policy_list,policy];
                end
            end
        end
    end
    Xcost(ii) = cost; %store the overall cost
    if ii<101 %store first 100 sample paths for visulization purposes
        path_x{ii}=xlist;
        path_y{ii}=ylist;
        d0_cell{ii}=d0_set;
        dmax_cell{ii}=dmax_set;
        policy_cell{ii} = policy_list;
    end
end
% t1 = toc
%% plotting the empirical CDF
prob_death = count_death/sample_size;

[f_1,x_1] = ecdf(Xcost,'Function','cdf','Alpha',0.05,'Bounds','on');
figure
plot(x_1,f_1,'b-','linewidth',2);
hold on
xline(median(Xcost),'b:','linewidth',2);
% xlim([2 7]);
ylim([0 1]);

ax = gca;
ax.FontSize = 10;
xlabel('overall cost (s)','Fontsize',15);
ylabel('probability of success','Fontsize',15);
lgd1 = legend('CDF','Median cost','location','northwest');
lgd1.FontSize = 18;
grid on;

%% sample path visualization
%-----------set up tenary coordiantes---------------
NN=800;
ratio = N/NN;
q=linspace(0,1,NN+1);
p=linspace(0,1,NN+1);
x1=[];x2=[];x3=[];

for i=1:1:NN+1
    x1ij=p(i)*ones(1,NN+1);
    x2ij=(1-p(i))*(1-q);
    x3ij=1-x1ij-x2ij;
    x1=[x1 x1ij];
    x2=[x2 x2ij];
    x3=[x3 x3ij];
 end

[x, y] = terncoords(x3, x2, x1);

K=NN;
Xtri=zeros(K+1,K+1); Ytri=zeros(K+1,K+1);
for j=1:K+1
    Xtri(j,1:K+1)= x(K*(j-1)+j:K*(j-1)+j+K);
    Ytri(j,1:K+1)= y(K*(j-1)+j:K*(j-1)+j+K);
end

%-------------select a sample path from stored ones-----------------------
indxx = 1;
ylist = path_y{indxx};
xlist = path_x{indxx};
d0_set = d0_cell{indxx};
dmax_set = dmax_cell{indxx};
policy_list = policy_cell{indxx};
this_cost = Xcost(indxx);

x1= ylist; x2=(1-ylist).*(1-xlist); x3=(1-ylist).*xlist;
[x, y] = terncoords(x3, x2, x1);

%find the indices where the treatment strategy switches from d=0 to d=dmax
%or vice versa
switchset = [1];
for i = 2:length(policy_list)
    if policy_list(i) ~= policy_list(i-1)
        switchset = [switchset i];
    end
end
switchset = [switchset length(policy_list)];

figure
hold on
% fill([0 1 0.5 0],[0 0 0.866 0],'w-','linewidth',1);
[AA,BB]=contourf(Xtri,Ytri,Dmat_det(1:ratio:end,1:ratio:end),[0,1]);
set(BB,'LineColor','none');

plot(linspace(0.4950,0.5050,3),[0.8574,0.8574,0.8574],'c:','linewidth',1.5);
plot(linspace(0.005,0.99,100),0.0087*ones(1,100),'c:','linewidth',1.5);

for i = 2:length(switchset)-1
    indxxx = switchset(i);
    if policy_list(indxxx) == 0
        plot(x(switchset(i-1):switchset(i)),y(switchset(i-1):switchset(i)),...
            'r-','markersize',2,'linewidth',1.5);
        
    else
        plot(x(switchset(i-1):switchset(i)),y(switchset(i-1):switchset(i)),...
            'g-','markersize',2,'linewidth',1.5);
        
    end
end

if policy_list(end) == 0
    plot(x(switchset(end-1):end),y(switchset(end-1):end),...
        'g-','markersize',2,'linewidth',1.5);
    
else
    plot(x(switchset(end-1):end),y(switchset(end-1):end),...
        'r-','markersize',2,'linewidth',1.5);
    
end

plot(x(1),y(1),'marker','o','MarkerFaceColor','m','MarkerEdgeColor','m','markersize',4.3);
hold off
axis equal
text(-0.09,0,'VOP','FontSize',11);
text(0.53,0.85,'GLY','FontSize',11);
text(1.01,0,'DEF','FontSize',11);
titlename = sprintf('The overall cost of this sample path is = %.3f',this_cost);
title(titlename);
set(gca,'visible','off');
set(findall(gca,'type','text'),'visible','on');

end




