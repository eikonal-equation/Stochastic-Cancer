function Xcost = CDF_threshold_aware_policy(Dmat_det,D_thres,Initial_budget,xloc,yloc,choice,sample_size)
%This function computes the CDF y = Pr(J <= s) that measures the
%probability of keeping the accumulative cost J under a given initial
%threshold/budget value s using thresold-aware policies (paper name)
%
%Dmat_det (input): the policy matrix of stationary policy computed from
%                  Gluzman et al.
%D_thres (input): stored threshold-aware policies on the entire (q,p,s)
%                 space (for every 0.005 cost)
%Initial_budget (input): the initial threshold/budget to work with
%xloc (input): x (or q)-coordinate of the starting tumor configuration
%yloc (input): y (or p)-coordinate of the starting tumor configuration
%choice (input): choice of the determination strategy, only 3 options:
%                (i) 'conservative'; (ii) 'aggressive'; (iii) 'majority'.
%sample_size(input): sample size of Monte Carlo simulations
%
%Xcost (output): Array of the accumulative (random) cost of each sample path
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

%% Monte Carlo initializations
if N >= 800
    % set dt adaptively so for every time step
    %(according to what we saved, we have policy data for every 0.005 budget)

    dt0 = 0.005/(s)/40; %travel 1/40-th of a slice when not using drugs
    dtmax = 0.005/(dmax+s)/10; %travel 1/10-th of a slice when at MTD rate

    %uniform discretization of thresholds on [0,Initial_budget]
    budget_list = 0:0.005:Initial_budget;

else % for demo purpose, if N < 800, we have policy data for every 0.02 budget

    dt0 = 0.02/(s)/40; %travel 1/40-th of a slice when not using drugs
    dtmax = 0.02/(dmax+s)/10; %travel 1/10-th of a slice when at MTD rate

    %uniform discretization of thresholds on [0,Initial_budget]
    budget_list = 0:0.02:Initial_budget;

end

M = length(budget_list); %number of s-slices

%create a global variable that has all the policies from s=0 to s=budget
d_interior = cell(1,M);
for i = 1:M
    d_interior{i} = D_thres(:,:,i)';
end

count_death=0; %counting number of deaths

%% Monte Carlo Main Loop
Xcost = zeros(1,sample_size); %initialize the cost-recording array for each sample
parfor ii = 1:sample_size
    cost=0; %accumulative cost
    dmax_set=[];
    d0_set=[1];
    policy_list = [];
    xlist=[xloc];
    ylist=[yloc];
    w1=[0];w2=[0];w3=[0];

    policy = 0; %initial policy, assuming we know it before testing
    policy_list=[policy_list,policy];
    dt=dt0; %the corresponding first time step

    %potential storage array of policies if running out of budget
    after_0 =[];
    after_max= [];

    %run EM for the first time step given the initial policy above
    j=2;
    w1(j)=w1(j-1)+sqrt(dt)*normrnd(0,1);
    w2(j)=w2(j-1)+sqrt(dt)*normrnd(0,1);
    w3(j)=w3(j-1)+sqrt(dt)*normrnd(0,1);
    %coordinates from last step
    xval=xlist(j-1);
    yval=ylist(j-1);
    %EM to approximate the location at the next step give the policy
    ylist(j) = yval+fp(xval,yval,policy)*dt+ yval*(1-yval)*s1*(w1(j)-w1(j-1))+...
        yval*(1-yval)*(1-xval)*s2*(w2(j)-w2(j-1))+...
        yval*(1-yval)*xval*s3*(w3(j)-w3(j-1));

    xlist(j) = xval+fq(xval,yval)*dt...
        + xval*(1-xval)*(s3*(w3(j)-w3(j-1))-s2*(w2(j)-w2(j-1)));

    cost = cost +(policy+s)*dt; %accumulating cost
    next_budget=Initial_budget-cost; %computing the remaining budget

    while ylist(j) > 0.01
        %find the indices of the current position on the grid
        kq=find(xlist(j)<=xx',1);
        kp=find(ylist(j)<=xx',1);

        %----------Determination of the policy at next step----------------
        if ((next_budget < 0) || abs(next_budget -0) < 10^(-5))
            % using the stationary policy computed from Gluzman et al.
            % after running out of budget
            if (kp == (ind_death)) %if we are just below the death barrier
                d3=Dmat_det(kp,kq-1);d4=Dmat_det(kp,kq);
                switch choice
                    case 'conservative'
                        %super-conservative determination strategy
                        if ((d3 + d4) == 2)
                            policy = dmax;

                            after_max=[after_max,j];
                            dt = dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'aggressive'
                        %super-aggressive determination strategy
                        if ((d3 + d4) > 0)
                            policy = dmax;

                            after_max=[after_max,j];
                            dt = dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'majority'
                        %vote by majority determination strategy
                        if ((d3 + d4) > 0)
                            policy = dmax;

                            after_max=[after_max,j];
                            dt = dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                end
            else
                %if in the interior
                d1=Dmat_det(kp-1,kq-1);d2=Dmat_det(kp-1,kq);
                d3=Dmat_det(kp,kq-1);d4=Dmat_det(kp,kq);
                d_square=[d1,d2,d3,d4];

                switch choice
                    case 'conservative'
                        %super-conservative determination strategy
                        if sum(d_square)== 4
                            policy = dmax;

                            after_max=[after_max,j];
                            dt = dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt = dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'aggressive'
                        %super-aggressive determination strategy
                        if sum(d_square) > 0
                            policy = dmax;

                            after_max=[after_max,j];
                            dt=dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt=dt0;
                        end
                        policy_list=[policy_list,policy];
                    case 'majority'
                        %vote by majority determination strategy
                        if sum(d_square)>=3
                            policy = dmax;

                            after_max=[after_max,j];
                            dt=dtmax;
                        else
                            policy = 0;

                            after_0=[after_0,j];
                            dt=dt0;
                        end
                        policy_list=[policy_list,policy];
                end
            end



        else % when we are still within the initial budget
            %finding the s-slices just above/below the remaining budget
            k = find(next_budget <= budget_list,1);
            d_bottom = d_interior{k-1};
            d_top    = d_interior{k};

            if kp == (ind_rec+1) %if just above the recovery barrier
                %we will only use the policy on the row that is just above
                %the recovery barrier
                d1=d_bottom(kp,kq-1);d2=d_bottom(kp,kq);
                d3=d_top(kp,kq-1);d4=d_top(kp,kq);
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

            else %if in the interior of the spatial domain
                %label the top corners as 1,2,3,4
                d1=d_bottom(kp-1,kq-1);d2=d_bottom(kp-1,kq);
                d3=d_bottom(kp,kq-1);d4=d_bottom(kp,kq);
                %label the bottom corners as 5,6,7,8
                d5=d_top(kp-1,kq-1);d6=d_top(kp-1,kq);
                d7=d_top(kp,kq-1);d8=d_top(kp,kq);
                d_cube=[d1,d2,d3,d4,d5,d6,d7,d8];

                switch choice
                    case 'conservative'
                        %super-conservative
                        if sum(d_cube)== 8
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
                        %super-aggressive
                        if sum(d_cube) > 0
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
                        %majority
                        if sum(d_cube) >= 6
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
        %------Run EM by applying the policy we just determined above------
        j=j+1;
        w1(j)=w1(j-1)+sqrt(dt)*normrnd(0,1);
        w2(j)=w2(j-1)+sqrt(dt)*normrnd(0,1);
        w3(j)=w3(j-1)+sqrt(dt)*normrnd(0,1);
        yval=ylist(j-1);xval=xlist(j-1);
        ylist(j) = yval+fp(xval,yval,policy)*dt+ yval*(1-yval)*s1*(w1(j)-w1(j-1))+...
            yval*(1-yval)*(1-xval)*s2*(w2(j)-w2(j-1))+...
            yval*(1-yval)*xval*s3*(w3(j)-w3(j-1));
        xlist(j) = xval+fq(xval,yval)*dt+xval*(1-xval)*(s3*(w3(j)-w3(j-1))-s2*(w2(j)-w2(j-1)));

        cost = cost+(policy+s)*dt; %accumulating cost

        next_budget = Initial_budget - cost; %computing the remaining budget


        %if we cross the failure barrier, stop
        if ylist(j) > 0.99
            count_death = count_death+1;

            %set the cost to be a large value (the theoretical value should be infinite)
            cost = 10^6;
            break
        end

    end
    Xcost(ii)=cost;
    %store the first 100 samples for visualization purpose
    if ii < 101
        path_x{ii} = xlist;
        path_y{ii} = ylist;
        d0_cell{ii} = d0_set;
        dmax_cell{ii}= dmax_set;
        after_0_cell{ii} = after_0;
        after_max_cell{ii} = after_max;
        policy_cell{ii} = policy_list;

    end
end


%% plotting the empirical CDF
prob_death = count_death/sample_size;

[f_1,x_1] = ecdf(Xcost,'Function','cdf','Alpha',0.05,'Bounds','on');
figure
plot(x_1,f_1,'b-','linewidth',2);
xlim([1 9]);
ylim([0 1]);

ax = gca;
ax.FontSize = 10;
xlabel('overall cost (s)','Fontsize',15);
ylabel('probability of success','Fontsize',15);
grid on;

%% sample path visualization
indxx = 1;
ylist = path_y{indxx};
xlist = path_x{indxx};
d0_set = d0_cell{indxx};
dmax_set = dmax_cell{indxx};
after_0 = after_0_cell{indxx};
after_max = after_max_cell{indxx};
policy_list = policy_cell{indxx};
this_cost = Xcost(indxx);

x1=ylist; x2=(1-ylist).*(1-xlist); x3=(1-ylist).*xlist;
[x, y] = terncoords(x3, x2, x1);

figure
hold on
fill([0 1 0.5 0],[0 0 0.866 0],'w-','linewidth',1);
plot(linspace(0.4950,0.5050,3),[0.8574,0.8574,0.8574],'c:','linewidth',1.5);
plot(linspace(0.005,0.99,100),0.0087*ones(1,100),'c:','linewidth',1.5);
plot(x(d0_set),y(d0_set),'g.','markersize',2.5,'linewidth',1.5);
plot(x(dmax_set),y(dmax_set),'r.','markersize',2.5,'linewidth',1.5);
if ~isempty(after_0)
    plot(x(after_0),y(after_0),'b.','markersize',2.5,'linewidth',1.5);
end
if ~isempty(after_max)
    plot(x(after_max),y(after_max),'y.','markersize',2.5,'linewidth',1.5);
end
plot(x(1),y(1),'marker','o','MarkerFaceColor','m','MarkerEdgeColor','m','markersize',4.3);
hold off
axis equal
text(-0.08,0,'VOP','FontSize',11);
text(0.53,0.85,'GLY','FontSize',11);
text(1.01,0,'DEF','FontSize',11);
titlename = sprintf('The overall cost of this sample path is = %.3f',this_cost);
title(titlename);
set(gca,'visible','off');
set(findall(gca,'type','text'),'visible','on');

end


