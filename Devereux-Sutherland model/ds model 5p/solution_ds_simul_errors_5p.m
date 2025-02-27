%--------------------------------------------------------------------------
% Two-country Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script computes Euler errors (up to fifth-order),
% replicating the results from Table 1 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: Euler errors')
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
approx = 3; %order of approximation

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------

%First, solve for the SSS
q1 = clock;
a_sss = compute_sss_ds_5p(model,params,M,eps_ind,0)
my_time1 = etime(clock, q1) 

%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
params(1) = a_sss;
%Deterministic steady-state
nxss=[ah1;af1;zeros(6,1)]; %for states
nyss = [1;1;z0;z0;1;1]; %for controls

%Check DSS is good
mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
%Compute derivatives
q2 = clock;
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);
my_time2 = etime(clock, q2)

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1; % evaluate at the model of interest (epsilon=1)
derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon is unit-root

T0 = 1000;
T = 10000;
%draw pseudo-random innovations
innovations = my_mvnrnd(zeros(n_e,1),Sigma,(T0 + (T-1)))'; % shocks from period 2 to T
%use Levintal's function simul.m to simulate the economy:
[yt,xt]=simul(x0,innovations,nyss,nxss,eta,derivs1,approx,0,model);
%[yt,xt]=simul_mod_pruning3(x0,innovations,nyss,nxss,eta,derivs1);


yt = yt(:,T0+1:end);
xt = xt(:,T0+1:end);

P = [betta gama kappa d0];
%Use monomials to discretize future innovations
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);
ln = length(weight_nodes);

my_errors =zeros(2,T);
%Loop to compute Euler errors
tic
for t=1:T
    
    my_errors(:,t) = log10(euler_errors_ds_5p(P,nxss,nyss,xt(:,t),epsi_nodes,weight_nodes,derivs1,eta,approx));
           
end
toc

errors_stats_cee1 = [mean(my_errors(1,:)) median(my_errors(1,:)) max(my_errors(1,:))]
errors_stats_pee = [mean(my_errors(2,:)) median(my_errors(2,:)) max(my_errors(2,:))]


