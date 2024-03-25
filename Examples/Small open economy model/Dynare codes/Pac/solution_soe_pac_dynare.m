%---------------------------------------------------------------------------
% Small open economy: Simulation calibrated with mexican data
% Auxiliary model is PAC
% This script computes Euler equation errors and kernel distributions
% Remember to run first the script pre_processing_soe_pac_full.m
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear
rng(0)

%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.2\matlab');
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

T0 = 1000;
T = 100000;
Sigma = M_.Sigma_e;
%draw pseudo-random innovations
innovations = mvnrnd([0 0],Sigma,(T0 + (T-1)))';

%DEFAULT: use Dynare matlab function simult_.m to simulate the economy, without pruning:
%myt =simult_(M_,options_,x0,mdr,innovations',3);
%
%ALTERNATIVE 1: same, but with pruning. You get nonsense because perturbation
%objects are treated as variables, so that the first-order component of the
%auxiliary model is a near unit-root process.
% options_.pruning = 1;
% myt =simult_(M_,options_,x0,mdr,innovations',3);
%
%ALTERNATIVE 2: Use simult_mod.m to simulate the
%economy, a slightly modified version of simult_ that implements pruning by
%treating the perturbation object sigma as a parameter.
options_.pruning = 1;
myt =simult_mod(M_,options_,x0,mdr,innovations',3);

myt = myt(:,T0+1:end);
xt = myt(1:4,:); %state variables
bt = xt(1,:); %bonds
ct = C0*myt(5,:); %consumption

%Compute and plot kernel density functions for NFA and consumption
pd_obj_b = fitdist(bt','Kernel','Kernel','epanechnikov');
XXb=-6:0.01:6;
%XX=0:0.01:70;
pd_b_points = pdf(pd_obj_b,XXb);
figure;plot(XXb,pd_b_points,'r');title('wealth');

pd_obj_c = fitdist(ct','Kernel','Kernel','epanechnikov');
XXc=0.85:0.01:1.15;
pd_c_points = pdf(pd_obj_c,XXc);
figure;plot(XXc,pd_c_points,'r');title('consumption');

%Compute and plot Euler errors
%inputs for the Euler errors function.
P2 = [betta*(1+psi2_calib) gama R0 C0];
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(2,Sigma);

lerrors =zeros(1,T);
tic
for t=1:T
    lerrors(t) =  log10(euler_errors_soe_dynare(P2,yss,mdr,myt(:,t),epsi_nodes,weight_nodes,3));
end
toc

errors_stats = ([mean(lerrors) median(lerrors) max(lerrors)])
figure;histogram(lerrors);title('errors CEE');

