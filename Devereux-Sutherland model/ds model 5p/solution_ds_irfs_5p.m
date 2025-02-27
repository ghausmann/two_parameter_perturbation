%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: Impulse responses.
% Standard version with asset holdings
% This script replicates Figure 3 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: impulse responses')
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
rho_d = 0.5; % persistence capital share 
uk_ul_corr = -0.5; %correlation between capital and labor income shocks 
z0 = betta*d0; %DSS asset price

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d];
%make eta
eta=[zeros(2,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];

%Moments
%Compute the cross moments:
n_e=4; %number of shocks.
M=gaussian_moments(n_e); %if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 0 uk_ul_corr 0;0 1 0 uk_ul_corr;uk_ul_corr 0 1 0;0 uk_ul_corr 0 1];
M.M2 = Sigma;

eps_ind = 8; %index of epsilon in state vector
%algo='gensylv';
algo='vectorize';
approx = 3; %order of approximation for dynamics

%color of plots?
my_color = 'b';
%my_color = 'r--';
%my_color = 'g:';

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------

%First, solve for the SSS
a_sss = compute_sss_ds_5p(model,params,M,eps_ind,0)
 
%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
params(1) = a_sss;

%Deterministic steady-state
nxss=[ah1;af1;zeros(6,1)]; %for states
nyss = [1;1;z0;z0;1;1]; %for controls

%Check DSS is good
mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
%Compute derivatives of policy rules
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1; % evaluate at the model of interest (epsilon=1)
derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon is unit-root
y0 = dr_gt(derivs1,nyss,approx,x0-nxss); %SSS of controls

%Shock to capital income share
x1 = x0;
x1(eps_ind-1) = 0.01;

T0 = 1;
T = 50;
Tr = 10;
time_h = 2;

%make innovations
innovations = zeros(4,(T0 + (T-1))); % shocks from period 2 to T

%use Levintal's function simul.m to simulate the economy:
[yt1,xt1]=simul(x1,innovations,nyss,nxss,eta,derivs1,3,0,model);
%uncomment this one for a second-order solution
%[yt1,xt1]=simul(x1,innovations,nyss,nxss,eta,derivs1,2,0,model);
%uncomment this one for a linear-corrected solution
%[yt1,xt1]=simul_linear_risk(x1,innovations,nyss,nxss,eta,derivs1,approx,1e-4);

yt = [y0 yt1];
xt = [x0 xt1];

%make variables of interest
zht = yt(3,1:end-1);
zft = yt(4,1:end-1);
aht = xt(1,2:end);
aft = xt(2,2:end);
wealth = zht.*aht + zft.*aft;

%Plot the results
time = 0:(Tr);

figure; 
subplot(2,2,1);
plot(time,zht(T0:T0+Tr),my_color,'LineWidth',2);
title('(a) Home asset price');
subplot(2,2,2);
plot(time,wealth(T0:T0+Tr),my_color,'LineWidth',2);
axis([0 Tr -0.5 0.5]);
title('(b) Wealth');
subplot(2,2,3);
plot(time,aht(T0:T0+Tr),my_color,'LineWidth',2);
title('(c) Home asset holdings');
subplot(2,2,4);
plot(time,aft(T0:T0+Tr),my_color,'LineWidth',2);
title('(c) Foreign asset holdings');


