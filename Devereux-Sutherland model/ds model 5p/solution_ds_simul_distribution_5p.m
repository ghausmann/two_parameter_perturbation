%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script replicates figure 4 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: ergodic distribution')
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
%------------------------------------------------------------------------
%kappa = 0; %Uncomment this instead to make the discount factor exogenous
%------------------------------------------------------------------------
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
%Compute derivatives
q2 = clock;
derivs0=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);
my_time2 = etime(clock, q2)
derivs1 = derivs0;

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1; % evaluate at the model of interest (epsilon=1)
derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon is unit-root

T0 = 1000;
T = 1000000;
%draw pseudo-random innovations
innovations = my_mvnrnd(zeros(n_e,1),Sigma,(T0 + (T-1)))'; % shocks from period 2 to T
%use Levintal's function simul.m to simulate the economy:
[yt1,xt1]=simul(x0,innovations,nyss,nxss,eta,derivs0,approx,0,model);
%use my modification to simulate with pruning
[yt2,xt2]=simul_mod_pruning3(x0,innovations,nyss,nxss,eta,derivs1);

yt1 = yt1(:,T0+1:end);
xt1 = xt1(:,T0+1:end);
zht1 = (yt1(3,1:end-1));
zft1 = (yt1(4,1:end-1));
aht1 = xt1(1,2:end);
aft1 = xt1(2,2:end);
wht1 = zht1.*aht1;
wft1 = zft1.*aft1;
Wealth1 = wht1 + wft1;
%Compute kernel density functions for NFA
pd_obj_b1 = fitdist(Wealth1','Kernel','Kernel','epanechnikov');
XXb1=min(Wealth1):0.01:max(Wealth1);
pd_b_points1 = pdf(pd_obj_b1,XXb1);

yt2 = yt2(:,T0+1:end);
xt2 = xt2(:,T0+1:end);
zht2 = (yt2(3,1:end-1));
zft2 = (yt2(4,1:end-1));
aht2 = xt2(1,2:end);
aft2 = xt2(2,2:end);
wht2 = zht2.*aht2;
wft2 = zft2.*aft2;
Wealth2 = wht2 + wft2;
pd_obj_b2 = fitdist(Wealth2','Kernel','Kernel','epanechnikov');
%XXb=-6:0.01:6;
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
    plot(XXb2,pd_b_points2,'r');
    title('(b) Exogenous discount factor');
end
    

