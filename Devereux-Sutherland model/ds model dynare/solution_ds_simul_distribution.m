%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script replicates figure 4 in the paper
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: ergodic distribution')
disp('-----------------------------------------------');

clear;

rng(0); %Fix seed for pseudorandom number generator.
%Add Dynare to the search path
addpath('C:\dynare\5.2\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

rho_y = 0.51;
M_.params(10) = rho_y;
kappa = 0.007; %Uzawa parameter
%kappa = 0; %Uzawa parameter
M_.params(9) = kappa;
zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------
%First, solve for the SSS
tic
a_sss = compute_sss_ds(M_,options_,oo_,eps_ind,0)
toc
% 
%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
M_.params(1) = a_sss;
yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

%Compute derivatives
tic
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);
toc

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
yss = mdr.ys;
y0 = yss; % start at the steady state
y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)
x0 = y0(5:12);
%SSS of all variables
y1 = dr_yt(mdr,yss,3,x0-yss(5:12,1),zeros(4,1));

T0 = 1000;
T = 1000000;
n_e=4; % number of shocks.
%draw pseudo-random innovations
innovations = mvnrnd(zeros(n_e,1),M_.Sigma_e,(T0 + (T-1)))'; % shocks from period 2 to T
%use Dynare's function simult_ to simulate the economy
myt1 = simult_(M_,options_,y1,mdr,innovations',3);
%use my modification of simult_ to simulate with pruning
myt2 = simult_mod(M_,options_,y1,mdr,innovations',3);

zht1 = myt1(1,1:end-1);
zft1 = myt1(2,1:end-1);
aht1 = myt1(5,1:end-1);
aft1 = myt1(6,1:end-1);
wht1 = zht1.*aht1;
wft1 = zft1.*aft1;
Wealth1 = wht1 + wft1;
%Compute kernel density functions for NFA
pd_obj_b1 = fitdist(Wealth1','Kernel','Kernel','epanechnikov');
XXb1=min(Wealth1):0.01:max(Wealth1);
pd_b_points1 = pdf(pd_obj_b1,XXb1);

zht2 = myt2(1,1:end-1);
zft2 = myt2(2,1:end-1);
aht2 = myt2(5,1:end-1);
aft2 = myt2(6,1:end-1);
wht2 = zht2.*aht2;
wft2 = zft2.*aft2;
Wealth2 = wht2 + wft2;
pd_obj_b2 = fitdist(Wealth2','Kernel','Kernel','epanechnikov');
XXb2=min(Wealth2):0.01:max(Wealth2);
pd_b_points2 = pdf(pd_obj_b2,XXb2);

%Plot the results
if kappa>0
    
    figure;
    subplot(1,2,1);
    plot(XXb1,pd_b_points1,'b');
    hold on;
    plot(XXb2,pd_b_points2,'r');
    title('(a) Endogenous discount factor');
    
    
elseif kappa==0
    subplot(1,2,2);
    plot(XXb1,pd_b_points1,'b');
    hold on;
    plot(XXb2,pd_b_points2,'r--');
    title('(b) Exogenous discount factor');
end
    

