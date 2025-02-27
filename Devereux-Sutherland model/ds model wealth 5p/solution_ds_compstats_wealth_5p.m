%--------------------------------------------------------------------------
% Two-country Devereux-Sutherland (2011) model: comparative statics.
% Version with real wealth, gross rates, and prices as states
% It replicates Figure 1, panel (a) of Appendix A.
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
gap = 0;

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa, rho_y, rho_d, gap];
%Moments
% Compute the cross moments:
n_e=4; % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 0 uk_ul_corr 0;0 1 0 uk_ul_corr;uk_ul_corr 0 1 0;0 uk_ul_corr 0 1];
M.M2 = Sigma;

eps_ind = 10; %index of epsilon in state vector
%algo='gensylv';
algo='vectorize';

%--------------------------------------------------------------------------
%Comparative statics: capital income share
%--------------------------------------------------------------------------

v_d = 0.1:0.025:1;
l_v_d = length(v_d);
sol_d = zeros(1,l_v_d);

q = clock;  

for t=1:l_v_d
    
    params_d = params;
    params_d(5) = v_d(t); % update std.
    % this function computes the SSS of nominal holdings of the home asset
    sol_2 = compute_sss_ds_wealth_5p(model,params_d,M,eps_ind,[0.2 betta*v_d(t)]);
    sol_d(t) = sol_2(1)/sol_2(2);
    
end

 my_time = etime(clock, q)

%Theoretical solution in DS paper
sol_ds_d = -0.5*(1 + ((1-v_d)./v_d)*uk_ul_corr*(tl/tk) );

%compare and plot the results
my_comp_d = [sol_d;sol_ds_d;sol_d-sol_ds_d];

figure;
subplot(1,2,1);
plot(v_d,sol_d(1,:),'LineWidth',2);
hold on;
plot(v_d,sol_ds_d,'r--','LineWidth',2);
title('(a) Comparative statics');
xlabel('\delta');
legend('SSS portfolio','Zero-order portfolio');

subplot(1,2,2);
title('(b) Dynamics')
xlabel('Time');

 