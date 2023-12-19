function Xcost = CDF_threshold_aware_policy_withDetOpt(Dmat_det,D_thres_filepath,Initial_budget,...
    Initial_policy,xloc,yloc,choice,sample_size,Total_Num_DataSlices,Sbar)
%This function computes the CDF y = Pr(J <= s) that measures the
%probability of keeping the accumulative cost J under a given initial
%threshold/budget value "s" using thresold-aware policies developed in
%the paper "Threshold-awareness in adaptive cancer therapy"
%(https://www.biorxiv.org/content/10.1101/2022.06.17.496649v2)
%(This one is for our 1st Example. 
% Here we use the deterministic policy after the budget runs out.)
%
%Dmat_det (input): the policy matrix of stationary policy computed from
%                  "Deterministic_Cancer_ValuePolicy_Ite.m"
%D_thres_filepath (input): file path for stored threshold-aware policies
%                          on the entire (q,p,s) space (for every 0.005 cost)
%Initial_budget (input): the initial threshold/budget to work with
%xloc (input): x (or q)-coordinate of the starting tumor configuration
%yloc (input): y (or p)-coordinate of the starting tumor configuration
%Initial_policy (input): the initial policy associated with the initial
%         configurations (assume to the known prior to simulations)
%choice (input): choice of the determination strategy, only 3 options:
%                (i) 'conservative'; (ii) 'aggressive'; (iii) 'majority'.
%sample_size(input): sample size of Monte Carlo simulations
%Total_Num_DataSlices(input): number of s-slices stored
%Sbar(input): captial Sbar implemented in the C++ code
%
%Xcost (output): Array of the accumulative (random) cost of each sample path
%
% Author: MingYi Wang, Cornell University
% Last modified: 12/2023
%
%% parameters
N = length(Dmat_det) - 1; %number of points along one side of spatial grid
dx = 1/N;
dy = dx;
% parameters from Gluzman et al.
% https://royalsocietypublishing.org/doi/10.1098/rspb.2019.2454
ba = 2.5;
bv = 2;
c = 1;
n = 4;

s1=0.15;s2=0.15;s3=0.15;%small diffusion constant
% s1=0.5;s2=0.5;s3=0.5;%large diffusion constant

s = 0.05;%treatment cost
dmax = 3;%MTD

xx = linspace(0,1,N+1);
r_b = 0.01; %recovery barrier
f_b = 1 - r_b; %death barrier
ind_rec = length(0:dx:r_b);%index for the stabilization barrier
ind_death = find(xx==f_b,1);%index for the death barrier


%drift functions (deterministic portion) of the cancer dynamics
fp = @(q,p,d) p*(1-p)*(ba/(n+1)-q*(bv-c)-d)-p*(s1^2*p-(s1^2*p^2+s2^2*(1-p)^2*(1-q)^2 ...
    + s3^2*(1-p)^2*q^2));
fq = @(q,p) q*(1-q)*(bv/(n+1)*(1+p+p^2+p^3+p^4)-c)+q*(1-q)*( (1-q)*s2^2-q*s3^2  );

D_thres = memmapfile(D_thres_filepath,'Format',{'uint8',[N+1,N+1,Total_Num_DataSlices],'dd'})
%% Monte Carlo initializations
% set dt adaptively so for every time step
ds = Sbar/(Total_Num_DataSlices-1);
dt0 = 0.005/(s)/40; %travel 1/40-th of a slice when not using drugs
dtmax = 0.005/(dmax+s)/10; %travel 1/10-th of a slice when at MTD rate

%uniform discretization of thresholds on [0,Initial_budget]
budget_list = 0:ds:Initial_budget;

count_death = 0; %counting number of deaths

%% Monte Carlo Main Loop
tic
Xcost = zeros(1,sample_size); %initialize the cost-recording array for each sample
parfor ii = 1:sample_size
    cost_set = [0]; %accumulative cost
    xlist = [xloc];
    ylist = [yloc];
    death = [];
    %Initialization of a standard 3D Brownian motion
    w1=[0];
    w2=[0];
    w3=[0];

    policy = Initial_policy; %initial policy, assuming we know it before testing
    policy_list = [policy];
    if (policy == 0)
        dmax_set = [];
        d0_set = [1];
        dt = dt0;
    else
        dmax_set = [1];
        d0_set = [];
        dt = dtmax;
    end

    %potential storage array of policies if running out of budget
    after_0 = [];
    after_max = [];

    j = 1;
    while ylist(j) >= r_b
        j = j+1; %update the time step
        %Sample Brownian motion as W_n - W_{n-1} = sqrt{dt}*N(0,1)
        w1(j) = w1(j-1)+sqrt(dt)*normrnd(0,1);
        w2(j) = w2(j-1)+sqrt(dt)*normrnd(0,1);
        w3(j) = w3(j-1)+sqrt(dt)*normrnd(0,1);

        %coordinates of the position from last step
        yval = ylist(j-1);
        xval = xlist(j-1);


        %EM to approximate the location at the next step give the policy
        ylist(j) = yval+fp(xval,yval,policy)*dt+ yval*(1-yval)*s1*(w1(j)-w1(j-1))+...
            yval*(1-yval)*(1-xval)*s2*(w2(j)-w2(j-1))+...
            yval*(1-yval)*xval*s3*(w3(j)-w3(j-1));

        xlist(j) = xval+fq(xval,yval)*dt...
            + xval*(1-xval)*(s3*(w3(j)-w3(j-1))-s2*(w2(j)-w2(j-1)));
        if xlist(j) > 1
            xlist(j) = 1;
        elseif xlist(j) < 0
            xlist(j) = 0;
        end

        %update the accumulative cost and remaining budget
        cost_set(j) = cost_set(j-1) + (policy+s)*dt; %accumulating cost
        next_budget = Initial_budget - cost_set(j); %find the remaining budget

        %if we cross the failure barrier, stop
        if ylist(j) > f_b
            count_death = count_death + 1;
            death = [death ii];
            %set the cost to be a large value (the theoretical value should be infinite)
            cost_set(j) = 10^6;
            break
        end
        %find the indices of the current position on the grid
        kq=find(xlist(j)<=xx',1);
        kp=find(ylist(j)<=xx',1);

        %----------Determination of the policy at next step----------------
        if ((next_budget < 0) || abs(next_budget -0) < 10^(-5))
            if (kp == (ind_death)) %if we are just below the death barrier
                d_2pt = [Dmat_det(kp,kq-1),Dmat_det(kp,kq)];
                policy_flag = Use_drug_or_not_2pt(d_2pt,choice);

            else
                %if in the interior
                d1=Dmat_det(kp-1,kq-1);d2=Dmat_det(kp-1,kq);
                d3=Dmat_det(kp,kq-1);d4=Dmat_det(kp,kq);
                d_square=[d1,d2,d3,d4];
                policy_flag = Use_drug_or_not_stat(d_square,choice);
            end

            if policy_flag
                policy = dmax;
                after_max=[after_max,j];
                dt = dtmax;

            else
                policy = 0;
                after_0=[after_0,j];
                dt = dt0;

            end
            policy_list = [policy_list,policy];

        else % when we are still within the initial budget
            %finding the s-slices just above/below the remaining budget
            k = find(next_budget <= budget_list,1);

            d1 = double(D_thres.Data.dd(kq-1,kp-1,k-1));
            d2 = double(D_thres.Data.dd(kq,kp-1,k-1));
            d3 = double(D_thres.Data.dd(kq-1,kp,k-1));
            d4 = double(D_thres.Data.dd(kq,kp,k-1));
            d5 = double(D_thres.Data.dd(kq-1,kp-1,k));
            d6 = double(D_thres.Data.dd(kq,kp-1,k));
            d7 = double(D_thres.Data.dd(kq-1,kp,k));
            d8 = double(D_thres.Data.dd(kq,kp,k));

            d_cube=[d1,d2,d3,d4,d5,d6,d7,d8];


            if kp == (ind_rec+1)
                %if just above the recovery barrier
                %we will only use the policy on the row that is just above
                %the recovery barrier
                d_square = [d3,d4,d7,d8];
                policy_flag = Use_drug_or_not_thres_square(d_square,choice);

            else

                policy_flag = Use_drug_or_not_thres(d_cube,choice);
            end

            if policy_flag
                policy = dmax;
                dmax_set=[dmax_set,j];
                dt = dtmax;

            else
                policy = 0;
                d0_set=[d0_set,j];
                dt = dt0;

            end
            policy_list = [policy_list,policy];

        end

    end
    Xcost(ii)= cost_set(end);
    %store the first 100 samples for visualization purpose
    if ii < 101
        path_x{ii} = xlist;
        path_y{ii} = ylist;
        d0_cell{ii} = d0_set;
        dmax_cell{ii}= dmax_set;
        after_max_cell{ii} = after_max;
        policy_cell{ii} = policy_list;
        cost_cell{ii} = cost_set;
    end
    
end
t_tol = toc

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
mycolor = [255/255, 10/255, 255/255];

indxx = 1;
y = path_y{indxx};
x = path_x{indxx};
d0_set = d0_cell{indxx};
dmax_set = dmax_cell{indxx};
after_max = after_max_cell{indxx};
after_0 = after_0_cell{indxx};
policy_list = policy_cell{indxx};
cost_list = cost_cell{indxx};
this_cost = Xcost(indxx);

after0_color = [255 153 51]/255;
aftermax_color = [102 51 0]/255;
after_start = [];

x1=y; x2=(1-y).*(1-x); x3=(1-y).*x;
[x, y] = terncoords(x3, x2, x1);

if ~isempty(after_0) || ~isempty(after_max)
    if isempty(after_0)
        after_start = after_max(1);
    elseif isempty(after_max)
        after_start = after_0(1);
    else

        after_start = min(after_0(1),after_max(1));
    end
    end_index1 = after_start;
    end_index2 = length(policy_list);
else
    end_index1 = length(policy_list);
end

switchset = [1];
for i = 2:end_index1
    if policy_list(i) ~= policy_list(i-1)
        switchset = [switchset i];
    end
end
switchset = [switchset end_index1];

if ~isempty(after_start)
    switch_after = after_start;
    for j = after_start+1:end_index2
        if policy_list(j) ~= policy_list(j-1)
            switch_after = [switch_after j];
        end
    end
    switch_after = [switch_after end_index2];
end

figure
hold on
fill([0 1 0.5 0],[0 0 0.866 0],'w-','linewidth',1);
for i = 2:length(switchset)-1
    indxxx = switchset(i);
    if policy_list(indxxx) == 0

        plot(x(switchset(i-1):switchset(i)),y(switchset(i-1):switchset(i)),...
            'r.-','markersize',2,'linewidth',1.5);
        hold on
    else

        plot(x(switchset(i-1):switchset(i)),y(switchset(i-1):switchset(i)),...
            'g.-','markersize',2,'linewidth',1.5);
        hold on
    end
end

if policy_list(end_index1) == 0
    plot(x(switchset(end-1):end_index1),y(switchset(end-1):end_index1),...
        'g.-','markersize',2,'linewidth',1.5);
    hold on
else
    plot(x(switchset(end-1):end_index1),y(switchset(end-1):end_index1),...
        'r.-','markersize',2,'linewidth',1.5);
    hold on
end

if ~isempty(after_start)
    for i = 2:length(switch_after)-1
        indxxx = switch_after(i);
        if policy_list(indxxx) == 0

            plot(x(switch_after(i-1):switch_after(i)),y(switch_after(i-1):switch_after(i)),...
                '.-','markersize',2,'linewidth',1.5,'Color',aftermax_color);
            hold on
        else

            plot(x(switch_after(i-1):switch_after(i)),y(switch_after(i-1):switch_after(i)),...
                '.-','markersize',2,'linewidth',1.5,'Color',after0_color);
            hold on
        end
    end

    if policy_list(end_index2) == 0
        plot(x(switch_after(end-1):end_index2),y(switch_after(end-1):end_index2),...
            '.-','markersize',2,'linewidth',1.5,'Color',after0_color);
        hold on
    else
        plot(x(switch_after(end-1):end_index2),y(switch_after(end-1):end_index2),...
            '.-','markersize',2,'linewidth',1.5,'Color',aftermax_color);
        hold on
    end
end
axis equal

plot(x(1),y(1),'marker','o','MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor,'markersize',4.5);
plot(linspace(0.4950,0.5050,3),[0.8574,0.8574,0.8574],'c:','linewidth',1.5);
plot(linspace(0.005,0.99,100),0.0087*ones(1,100),'c:','linewidth',1.5);

titlename = sprintf('The overall cost of this sample path is = %.3f',this_cost);
title(titlename);
end

%% subfunctions
function policy_flag = Use_drug_or_not_stat(neighbors,choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 4
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'majority'
        if sum(neighbors) >= 2
            policy_flag = true;
        else
            policy_flag = false;
        end
end
end


function policy_flag = Use_drug_or_not_thres(neighbors,choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 8
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'majority'
        if sum(neighbors) >= 5
            policy_flag = true;
        else
            policy_flag = false;
        end

end
end

function policy_flag = Use_drug_or_not_thres_square(neighbors,choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 4
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'majority'
        if sum(neighbors) >= 3
            policy_flag = true;
        else
            policy_flag = false;
        end

end
end


function policy_flag = Use_drug_or_not_2pt(neighbors,choice)
switch choice
    case 'conservative'
        if sum(neighbors) == 2
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'aggressive'
        if sum(neighbors) > 0
            policy_flag = true;
        else
            policy_flag = false;
        end
    case 'majority'
        if sum(neighbors) > 0
            policy_flag = true;
        else
            policy_flag = false;
        end

end
end




