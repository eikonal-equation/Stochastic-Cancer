function Xcost = CDF_stationary_policy_EG2(Dmat_det,xloc,yloc,Initial_policy,choice,sample_size)
%This function computes the CDF y = Pr(J <= s) that measures the
%probability of success where the cumulative cost J is within any positive
%threshold value "s" using the stationary policy computed with the base
%model from Carrère 2017 https://doi.org/10.1016/j.jtbi.2016.11.009
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
%
% Author: MingYi Wang, Cornell University
% Last modified: 12/2023
%
%% parameters
N = length(Dmat_det)-1; %number of points along one side of spatial grid
dx = 1/N;

%parameters from Carrere 2017 https://doi.org/10.1016/j.jtbi.2016.11.009
time_factor = 1;
gs = 0.031*time_factor;
gr = 0.026*time_factor;
m = 30;
K = 0.48; %of order 1e7
beta = 6.25*time_factor; %of order 1e-7
beta_hat = K*beta;
alpha = 0.06*time_factor;

s1 = 0.15*sqrt(time_factor);
s2 = s1;
% s1=0.5;s2=0.5;%larger diffusion constant

s=0.05;%treatment cost
dmax=3;%MTD

xx = linspace(0,1,N+1);
r_b = 0.01; %recovery barrier
f_b = 1 - r_b; %death barrier
ind_rec = length(0:dx:r_b);%index for the recovery barrier
ind_death = find(xx==f_b,1);%index for the death barrier

%drift functions (deterministic portion) of the cancer dynamics
fp = @(q,p,d) p*(1-p)*(gs*q+gr*(1-q)) - alpha*q*p*d - beta_hat*p^2*q*(1-q);
fq = @(q,p,d) (1-p)*(1-q)*q*(gs-gr) - alpha*q*(1-q)*d + beta_hat*p*q^2*(1-q)...
    +(1-p)^2*q*(1-q)*(s2^2*(1-q) - s1^2*q + s1*s2);

%diffusion functions (stochastic portion) of the cancer dynamics
fpp = @(q,p) p*(1-p)*(s1*q+s2*(1-q));
fqq = @(q,p) (1-p)*q*(1-q)*(s1-s2);

%% Monte Carlo Simulations
%set the time step adaptively such that for each step, we spend
% dt0 = 0.01/(s)/40; %1/40-th of 0.01 cost when not using drugs
% dtmax = 0.01/(dmax+s)/10; %1/10-th of 0.01 cost when using drugs at MTD rate
count_death = 0; %initialize the number of deaths to 0
dt0 = 0.001;
dtmax = 0.001;
% tic
parfor ii = 1:sample_size
    cost_set = [0]; %accumulative cost
    ylist = [yloc];
    xlist = [xloc];
    % 1D Brownian motion
    w1 = [0];

    
    %initial policy (assuming known)
    policy = Initial_policy;
    if (policy == 0)
        dmax_set = [];
        d0_set = [1];
        dt = dt0;
    else
        dmax_set = [1];
        d0_set = [];
        dt = dtmax;
    end
    policy_list = [policy];

    j = 1;
    while ylist(j) >= r_b
        j = j+1; %update the time step
        %Sample Brownian motion as W_n - W_{n-1} = sqrt{dt}*N(0,1)
        w1(j) = w1(j-1)+sqrt(dt)*normrnd(0,1);

        %coordinates of the position from last step
        yval = ylist(j-1);
        xval = xlist(j-1);
        %diffusion terms
        x_diff = fqq(xval,yval);
        y_diff = fpp(xval,yval);
        %EM to approximate the location at the next step give the policy
        ylist(j) = yval + fp(xval,yval,policy)*dt + y_diff*(w1(j)-w1(j-1));

        xlist(j) = xval + fq(xval,yval,policy)*dt + x_diff*(w1(j)-w1(j-1));

        if xlist(j) > 1
            xlist(j) = 1;
        elseif xlist(j) < 0
            xlist(j) = 0;
        end

        cost_set(j) = cost_set(j-1) + (policy+s)*dt; %accumulating cost

        %if we cross the failure barrier, stop
        if ylist(j) > f_b
            count_death = count_death + 1;
            %set the cost to be a large value (the theoretical value should be infinite)
            cost_set(j) = 10^6;
            break
        end

        if ylist(j) >= r_b
            %find the indices of the current position on the grid
            kq=find(xlist(j)<=xx',1);
            kp=find(ylist(j)<=xx',1);

            %----------Determination of the policy at next step----------------
            if (kp == (ind_death)) %if we are just below the death barrier
                d3 = Dmat_det(kp,kq-1);
                d4 = Dmat_det(kp,kq);
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
    Xcost(ii) = cost_set(end); %store the overall cost
    if ii<101 %store first 100 sample paths for visulization purposes
        path_x{ii}=xlist;
        path_y{ii}=ylist;
        d0_cell{ii}=d0_set;
        dmax_cell{ii}=dmax_set;
        policy_cell{ii} = policy_list;
        cost_cell{ii} = cost_set;
    end
    
end
%% plotting the empirical CDF
prob_death = count_death/sample_size;

[f_1,x_1] = ecdf(Xcost,'Function','cdf','Alpha',0.05,'Bounds','on');
figure
plot(x_1,f_1,'b-','linewidth',2);
hold on
xline(median(Xcost),'b:','linewidth',2);
xlim([30 120]);
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
if N >= 800 %if N >= 800, we suggest downsampling it to 801x801
    NN = 800;
    q=linspace(0,1,NN+1);
    p=linspace(0,1,NN+1);

else
    q=linspace(0,1,N+1);
    p=linspace(0,1,N+1);
end
[X,Y] = meshgrid(q,p);
%-------------select a sample path from stored ones-----------------------
N_plot = length(X) - 1;
indxx = 1;
y = path_y{indxx};
x = path_x{indxx};
d0_set = d0_cell{indxx};
dmax_set = dmax_cell{indxx};
policy_list = policy_cell{indxx};
this_cost = Xcost(indxx);

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
[AA,BB]=contourf(X,Y,Dmat_det(1:N/N_plot:end,1:N/N_plot:end),[0 1]);
set(BB,'LineColor','none');
yline(r_b,'c:','LineWidth',2);
yline(f_b,'c:','LineWidth',2);

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
plot(x(1),y(1),'marker','o','MarkerFaceColor','m','MarkerEdgeColor','m','markersize',5);
hold off
xlabel('Fraction of the Sensitive (Q)','FontSize',14);
ylabel('Total size (P)','FontSize',14);
axis equal
titlename = sprintf('The overall cost of this sample path is = %.3f',this_cost);
title(titlename);

end





