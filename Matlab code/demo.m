%% Optimal cancer therapy under stochasticity Demo file
%
% This script is usued to demonstrate how to parse date from running our
% C++ code into Matlab and run our Matlab code to generate CDFs and the
% deterministic-optimal policy.
%
% Due to demo purpose only, we will only use a coarse grid. Thus, the
% results will NOT achieve high enough accuracy.
% If you want to reproduce the results we show in our paper, please see
% section 5S of our Supplementary Materials for implementation details.
%
% Author: MingYi Wang, Cornell University
% Last modified: 10/24/22
%
clear all;
close all;
%% Generating the deterministic-optimal policy

% we've used a 3201x3201 uniform grid on the unit qp-square in our actual
% implementation. Here we are using a 401x401 grid as demo
N_det = 400;

% Call our given function to generate the deterministic-optimal policy
tic
[U,Dmat_det] = Deterministic_Cancer_ValuePolicy_Ite(N_det);
t1 = toc;
fprintf('The elapsed time of generating the det-optimal policy is %.3f seconds.\n',t1)
%% Parse data file from running C++ code

% First, open the data file. We only provide a 401x401x301 threshold-aware
% policy data file for demonstration.
% In our actual implementation, we've always generated a 1601x1601x1201
% policy date file.
tic
fileId = fopen('test_policy.dat');
D_MultArr = fread(fileId,'uint8'); % we've stored the optimal policy as either 0 or 1 ("bool" in C++)
D_MultArr = reshape(D_MultArr,[401,401,301]); % reshape it to be a 401x401x301 multi-dimensional array
t2 = toc;
fprintf('The elapsed time of loading the policy data file is %.3f seconds.\n',t2)

%% Generating CDFs with the deterministic-optimal policy
% picking a starting location such that it is initially in the blue
% (drug-off) region
xloc = 0.26;
yloc = 0.665;
% picking the policy determination strategy
choice = 'conservative';
% setting the sample size for MC simulation. We will use 1000 samples for
% demo purpose only. In our actual implementation, we've used 10^5 samples.
sample_size = 10^3;
% Call our function to generate the CDF
tic
Xcost_stationary = CDF_stationary_policy(Dmat_det,xloc,yloc,choice,sample_size);
t4 = toc;
fprintf('The elapsed time of generating the CDF with det-optimal policy is %.3f seconds.\n',t4)

%% Generating CDFs with threshold-aware policies
% Using the same initial tumor configuration, the same policy determination
% startegy, and the same sample size as above.
%
% Set up the initial budget sbar
sbar = 5;
% Call our function to generate the CDF
tic
Xcost_thres = CDF_threshold_aware_policy(Dmat_det,D_MultArr,sbar,xloc,yloc,choice,sample_size);
t5 = toc;
fprintf('The elapsed time of generating the CDF with threshold-aware policy is %.3f seconds.\n',t5)
