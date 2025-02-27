%---------------------------------------------------------------------------
% Small open economy: Simulation calibrated with mexican data
% Auxiliary model is PAC
% This script performs a stochastic simulation of the model.
% Compare the output with what you get using "solution_1p_soe.m" in the folder "one_param_soe". 
% Remember to run first the script pre_processing_soe_pac_full.m
%--------------------------------------------------------------------------

clear
rng(0); %fix seed for full replication

%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
%Load the pre-processing data
load('my_soe_pac_full.mat');

%index of the epsilon variable in the vector yt (Dynare's endogenous
%variables)
eps_ind = 4;

%--------------------------------------------------------------------------
% Calibration and perturbation solution
%--------------------------------------------------------------------------

%Calibrate pssi2, and indirectly the discount factor of the model of
%interest:
psi2_calib = compute_betta_soe_pac_dynare(b_bar,M_,options_,oo_,eps_ind,0)
M_.params(12) = psi2_calib;
%Just a check
b1_check = compute_sss_soe_pac_dynare(M_,options_,oo_,eps_ind,0)

%Compute the perturbation solution:
%The external algorithm is Dynare matlab function resol.m
[mdr, info, M_, oo_] = resol(0, M_, options_, oo_);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------

yss = mdr.ys;
x0 = yss; % start at the steady state
%Impose the model of interest (epsilon=1)
x0(eps_ind)=1;

T0 = 100;
T = 1000;
Sigma = M_.Sigma_e;
%draw pseudo-random innovations
innovations = my_mvnrnd([0 0],Sigma,(T0 + (T-1)))';

%DEFAULT: use Dynare matlab function simult_.m to simulate the economy, without pruning:
myt =simult_(M_,options_,x0,mdr,innovations',3);
myt = myt(:,T0+1:end);
xt = myt(1:4,:); %state variables
bt = xt(1,:); %bonds
ct = C0*myt(5,:); %consumption

figure;plot(bt);

