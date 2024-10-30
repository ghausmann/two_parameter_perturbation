%---------------------------------------------------------------------------
% One bond, two-country model: Comparative statics and simulation
% This script replicates Figure 10 in Appendix D.4
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%---------------------------------------------------------------------------

clear
rng(0); %fix seed for full replication
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
betta = 0.9851; %discount factor of the auxiliary model
gama = 2; %risk aversion
rho_y = 0.749; %auto-corr income
rho_eps = 1;
ty = 0.02723*((1-rho_y^2)^0.5); %conditional std. of income shocks
b_bar = 0; %Target for sss of NFA
gap = 0; %just a initial value (to be calibrated)
psi1 = 0.00001; %improves accuracy

C0h = 1 + (1-betta)*b_bar;
C0f = 1 - (1-betta)*b_bar;

% Make a vector of parameter values:
params = [betta, gama, rho_y, rho_eps, ty, b_bar, C0h, C0f, gap, psi1];
% Make eta
eta=[0 0;ty 0;0 ty;0 0];

%Moments
% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
Sigma = [1 0;0 1];
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(n_e,Sigma);
% Choose algorithm for the Sylvester equation
algo='gensylv'; % Kamenik algorithm
eps_ind = 4; %index of epsilon in state vector
approx_s = 2; %order of approximation for SSS
approx_d = 3; %order of approximation for dynamics

%--------------------------------------------------------------------------
% Comparative statics: 
%--------------------------------------------------------------------------
m_gap = 0:0.00001:0.0001;
l_gap = length(m_gap);

bsss = zeros(l_gap,1);
Rsss = zeros(l_gap,1);
sss_cond = zeros(l_gap,1);
Chsss = zeros(l_gap,1);
Cfsss = zeros(l_gap,1);
v_y1 = zeros(3,l_gap);
guessb = 0;

for t=1:l_gap
    
    params_cs = params;
    params_cs(9) = m_gap(t);
    % this function computes the SSS of bonds
    bsss(t) = compute_sss_bond2(model,params_cs,M,eta,eps_ind,approx_s,guessb);
    guessb = bsss(t);
    b_bar1 = bsss(t);
    
    C0h1 = 1 + (1-betta)*b_bar1;
    C0f1 = 1 - (1-betta)*b_bar1;
    
    params_cs(6) = b_bar1;
    params_cs(7) = C0h1;
    params_cs(8) = C0f1;
    
    nxss1 = [b_bar1;0;0;0];
    nyss1 = [1/betta;1;1];
    
    %Compute derivatives of policy rules
    derivs1=solve_dsge(model,params_cs,M,eta,nxss1,nyss1,approx_s,algo);
    
    %Compute the SSS of consumptions
    x0 = nxss1;
    x0(eps_ind) = 1; % impose the model of interest
    derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon follows a unit-root process.
    y1=dr_gt(derivs1,nyss1,approx_s,(x0-nxss1));
    v_y1(:,t) = y1;
    Rsss(t) = y1(1);
    sss_cond(t) = (betta+m_gap(t))*Rsss(t);
    Chsss(t) = C0h1*y1(2);
    Cfsss(t) = C0f1*y1(3);
    
end

%--------------------------------------------------------------------------
% Dynamics
%--------------------------------------------------------------------------
%Calibrate the gap between countries' impatience consistent with a SSS NFA
%target for the Home country.
b_target = 1;
mysol_calib = compute_calib_bond2(model,params,M,eta,eps_ind,approx_s,0,b_target)

b_bar1 = b_target;
gap = mysol_calib;
params(9) = gap;

C0h = 1 + (1-betta)*b_bar1;
C0f = 1 - (1-betta)*b_bar1;

params(6) = b_bar1;
params(7) = C0h;
params(8) = C0f;

nxss = [b_bar1;0;0;0];
nyss = [1/betta;1;1];

%Compute the perturbation solution:
%First, verify the DSS
dss_check=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
%The external algorithm is Levintal's function solve_dsge.m
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx_d,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------

x0=nxss; % start at the steady state
%Impose the model of interest (epsilon=1)
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1;

%Make sure the solution is theoretically valid:
y1=dr_gt(derivs,nyss,approx_d,(x0-nxss))
sss_cond = y1(1)*(betta+gap)
sss_cond_diff = 1 - sss_cond

T0 = 1000;
T = 1000000;
Te = 100000;
%draw pseudo-random innovations
innovations = mvnrnd([0 0],Sigma,(T0 + (T-1)))';

%DEFAULT: Use Levintal's function simul.m to simulate the economy, without pruning:
%[yt,xt]=simul(x0,innovations,nyss,nxss,eta,derivs,approx_d,0,model);
%
% ALTERNATIVE 1: same, but with pruning. You will get inaccurate results
% because perturbation objects are treated as variables, so that the
% first-order component of the auxiliary model is a near unit-root process.
[yt,xt]=simul(x0,innovations,nyss,nxss,eta,derivs,approx_d,1,model);
%
%ALTERNATIVE 2: Use my own function simul_mod_pruning3.m to simulate the
%economy, a heavily modified version of simul.m that implements pruning by
%treating perturbation objects as parameters (only for third-order).
%[yt,xt]=simul_mod_pruning3(x0,innovations,nyss,nxss,eta,derivs);

yt = yt(:,T0+1:end);
xt = xt(:,T0+1:end);

bht = xt(1,:);
yht = exp(xt(2,:));
yft = exp(xt(2,:));
cft = C0f*yt(3,:);

%Compute and plot kernel density functions for NFA and consumption
pd_obj_b = fitdist(bht(1,:)','Kernel','Kernel','epanechnikov');
XXb=min(bht):0.01:max(bht);
%XX=0:0.01:70;
pd_b_points = pdf(pd_obj_b,XXb);

%---------------------------------------------------------------------------
%Compute Euler equation errors:
%---------------------------------------------------------------------------

%inputs for the Euler errors function.
P = [betta;gama;C0h;C0f;gap];
lerrors =zeros(2,Te);
for t=1:Te

    lerrors(:,t) =  log10(euler_errors_bond2(P,nxss,nyss,xt(:,t),epsi_nodes,weight_nodes,derivs,eta,approx_d));

end

errors_stats = ([mean(lerrors(1,:)) median(lerrors(1,:)) max(lerrors(1,:));mean(lerrors(2,:)) median(lerrors(2,:)) max(lerrors(2,:))])

figure;
subplot(1,2,1);
plot(m_gap,Chsss,'b');
hold on;
plot(m_gap,Cfsss,'r--');
title('SSS consumptions');

subplot(1,2,2);
plot(XXb,pd_b_points);
title('Distribution of Home NFA');


