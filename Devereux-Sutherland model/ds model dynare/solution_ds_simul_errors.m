%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script computes Euler errors (up to third-order)
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: Euler errors')
disp('-----------------------------------------------');

clear;

%Add Dynare to the search path
addpath('C:\dynare\5.2\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

rho_y = 0.51;
M_.params(10) = rho_y;
zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------
%First, solve for the SSS

tic;a_sss = compute_sss_ds(M_,options_,oo_,eps_ind,0);toc

% 
%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
M_.params(1) = a_sss;
%Deterministic steady-state
yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

%Compute derivatives

tic;[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);toc


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
T = 10000;
n_e=4; % number of shocks.
%draw pseudo-random innovations
innovations = mvnrnd(zeros(n_e,1),M_.Sigma_e,(T0 + (T-1)))'; % shocks from period 2 to T
%use Dynare's function simult_ to simulate the economy
myt =simult_(M_,options_,y1,mdr,innovations',3);
xt = myt(5:12,:); %states
yt = [myt(1:4,:);myt(13:end,:)]; %controls

P = [betta gama kappa d0];
%Use monomials to discretize future innovations
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(4,M_.Sigma_e); %monomials to approximate expectations
%Loop to compute Euler errors
my_errors =zeros(2,T);

tic
for t=1:T
    
    my_errors(:,t) = log10(euler_errors_ds(P,yss,mdr,myt(:,t),epsi_nodes,weight_nodes,3));

end
toc

errors_stats_cee1 = [mean(my_errors(1,:)) median(my_errors(1,:)) max(my_errors(1,:))]
errors_stats_pee = [mean(my_errors(2,:)) median(my_errors(2,:)) max(my_errors(2,:))]

