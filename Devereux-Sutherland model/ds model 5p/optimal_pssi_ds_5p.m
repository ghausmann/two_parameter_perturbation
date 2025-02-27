%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: optimization routine to set
% pssi (free parameter that controls PAC-like modifications)
% Standard version with asset holdings
%
% This routine is discussed in Section 4.2 of the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: simulation')
disp('-----------------------------------------------');

clear;

rng(0); %Fix seed for pseudorandom number generator.
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
a = 0; %PAC location parameter
betta = 0.96; %discount factor
gama = 2; %risk-aversion
pssi = 5e-6; %PAC control parameter
d0 = 0.167; %Share of capital income
tk = 0.018; %std. capital income
tl = 0.018; %std. labor income
rho_eps = 1; % persistence of epsilon
kappa = 0.007; %Uzawa parameter
rho_y = 0.51; % persistence income
rho_d = 0; % persistence capital share 
uk_ul_corr = -0.5; %correlation between capital and labor income shocks 
z0 = betta*d0; %DSS asset price

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d];
%make eta
eta=[zeros(2,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];

%Moments
%Compute the cross moments:
n_e=4; % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 0 uk_ul_corr 0;0 1 0 uk_ul_corr;uk_ul_corr 0 1 0;0 uk_ul_corr 0 1];
M.M2 = Sigma;

eps_ind = 8; %index of epsilon in state vector
%algo='gensylv';
algo='vectorize';
approx = 3;

%--------------------------------------------------------------------------
% Optimization routine
%--------------------------------------------------------------------------

[~,epsi_nodes,weight_nodes] = Monomials_2(4,Sigma);
P = [betta gama kappa d0 z0];

%Objective function
myfun=@(x)errors_routine_5p(x,model,params,eta,M,eps_ind,approx,epsi_nodes,weight_nodes,P);
options = optimoptions('lsqnonlin','Display','iter');
t =clock; 
%Use the built-in function lsqnonlin to minimize
my_est0 = lsqnonlin(myfun,log(1e-5),[],[],options);
pssi_star = exp(my_est0)
my_time1 = etime(clock,t)
