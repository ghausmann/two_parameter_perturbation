%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: comparative statics.
% Standard version with asset holdings
% This script replicates Figure 2 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: comparative statics')
disp('-----------------------------------------------');

clear;

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
rho_y = 0; % persistence income
rho_d = 0; % persistence capital share 
uk_ul_corr = -0.5; %correlation between capital and labor income shocks 
z0 = betta*d0; %DSS asset price

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d];

%Moments
%Compute the cross moments:
n_e=4; %number of shocks.
M=gaussian_moments(n_e); %if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 0 uk_ul_corr 0;0 1 0 uk_ul_corr;uk_ul_corr 0 1 0;0 uk_ul_corr 0 1];
M.M2 = Sigma;

eps_ind = 8; %index of epsilon in state vector
%Choose algorithm for first-order solution
%algo='gensylv';
algo='vectorize';
approx = 2; %order of approximation for dynamics

%--------------------------------------------------------------------------
%Comparative statics: capital income share
%--------------------------------------------------------------------------
v_d = 0.1:0.025:1;
l_v_d = length(v_d);
sol_d = zeros(1,l_v_d);

for t=1:l_v_d
    
    params_d = params;
    params_d(5) = v_d(t); % update income share
    % this function computes the SSS of nominal holdings of the home asset
    sol_d(t) = compute_sss_ds_5p(model,params_d,M,eps_ind,0);
        
end

%Theoretical solution in DS paper
sol_ds_d = -0.5*(1 + ((1-v_d)./v_d)*uk_ul_corr*(tl/tk) );
%compare results
my_comp_d = [sol_d;sol_ds_d;sol_d-sol_ds_d];

%--------------------------------------------------------------------------
%Comparative statics: correlation between capital and income shocks
%--------------------------------------------------------------------------
v_c = -1:0.1:1;
l_v_c = length(v_c);
sol_c = zeros(1,l_v_c);

for t=1:l_v_c
    
    uk_ul_corr_c = v_c(t);
    Sigma_c = [1 0 uk_ul_corr_c 0;0 1 0 uk_ul_corr_c;uk_ul_corr_c 0 1 0;0 uk_ul_corr_c 0 1];
    Mc = M;
    Mc.M2 = Sigma_c;
    sol_c(t) = compute_sss_ds_5p(model,params,Mc,eps_ind,0);
        
end

%Theoretical solution in DS paper
sol_ds_c = -0.5*(1 + ((1-d0)/d0)*v_c*(tl/tk) );
%compare results
my_comp_c = [sol_c;sol_ds_c;sol_c-sol_ds_c];

%Make figures
figure;
subplot(1,2,1);
plot(v_d,sol_d,'LineWidth',2);
hold on;
plot(v_d,sol_ds_d,'r--','LineWidth',2);
title('(a) As a function of $\overline{\delta}$','Interpreter','latex');
xlabel('$\overline{\delta}$','Interpreter','latex');
legend('SSS portfolio','Zero-order portfolio');

subplot(1,2,2);
plot(v_c,sol_c,'LineWidth',2);
hold on;
plot(v_c,sol_ds_c,'r--','LineWidth',2);
title('(b) As a function of $\mu_{KL}$','Interpreter','latex');
xlabel('$\mu_{KL}$','Interpreter','latex');
legend('SSS portfolio','Zero-order portfolio');

