%--------------------------------------------------------------------------
% Two-country Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script replicates Figure 5 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: SSS matters')
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
d0 = 1; %Share of capital income (NOTE: only capital income shocks active)
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
approx = 3; %approximation order

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------

%First, solve for the SSS
a_sss = compute_sss_ds_5p(model,params,M,eps_ind,0)

%Vector of approximation points for asset holdings
a_v = [-0.8:0.1:-0.6 -0.59:0.01:-0.41 -0.4:0.1:-0.20];
l_a_v = length(a_v);

%Vectors of moments (to be filled in the coming loop)
errors_cee_v = zeros(l_a_v,1);
errors_pee_v = zeros(l_a_v,1);
my_corr = zeros(l_a_v,1);
my_std = zeros(l_a_v,1);

%draw pseudo-random innovations
T0 = 1;
T = 1000;
innovations = my_mvnrnd(zeros(n_e,1),Sigma,(T0 + (T-1)))'; % shocks from period 2 to T
%Compute Euler errors for the given simulation
P = [betta gama kappa d0 tk tl];
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);
sigma_nodes = [tk tk tl tl].*epsi_nodes;
ln = length(weight_nodes);


for n = 1:l_a_v
    
    a_sss = a_v(n); %fix the approximation point
    %Recalculate DSS of auxiliary model
    ah1 = a_sss;
    af1 = -a_sss;
    params(1) = a_sss;
        
    %Deterministic steady-state
    nxss=[ah1;af1;zeros(6,1)]; %for states
    nyss = [1;1;z0;z0;1;1]; %for controls
    
    %mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
    derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);
        
    %--------------------------------------------------------------------------
    % Simulation of the model of interest
    %--------------------------------------------------------------------------
    x0=nxss; % start at the steady state
    x0(eps_ind) = 1;
    derivs.hx(eps_ind,eps_ind) = 1; % evaluate at the model of interest (epsilon=1)
    
    %use Levintal's function simul.m to simulate the economy:
    [yt0,xt0]=simul(x0,innovations,nyss,nxss,eta,derivs,approx,0,model);
    %[yt0,xt0]=simul_mod_pruning3(x0,innovations,nyss,nxss,eta,derivs);
       
    yt = yt0(:,T0+1:end);
    xt = xt0(:,T0+1:end);
    ct = yt(1,:);
    cft = yt(2,:);
    
    %Compute second moments
    my_smoments = corrcoef(ct,cft);
    my_corr(n) = my_smoments(1,2);
    my_std(n) = std(xt(1,:));
    
    %Loop to compute Euler errors    
    my_errors =zeros(2,T);
    for t=1:T
        my_errors(:,t) = (euler_errors_ds_5p(P,nxss,nyss,xt(:,t),epsi_nodes,weight_nodes,derivs,eta,approx));
    end
    
    %Compute mean Euler errors
    errors_cee_v(n) = (mean(my_errors(1,:)));
    errors_pee_v(n) = (mean(my_errors(2,:)));
       
end

%Plot the results
figure;
subplot(1,2,1);
plot(a_v,errors_cee_v);
hold on;
plot(a_v,errors_pee_v','r--');
xlabel('$\bar{a}$','Interpreter','latex');
title('Euler equation errors');
legend('Consumption eq.','Portfolio eq.')
subplot(1,2,2);
plot(a_v,my_corr);
hold on;
plot(a_v,my_std','r--');
xlabel('$\bar{a}$','Interpreter','latex');
title('Second moments');
legend('Corr(ct,cst)','std(aht)')


