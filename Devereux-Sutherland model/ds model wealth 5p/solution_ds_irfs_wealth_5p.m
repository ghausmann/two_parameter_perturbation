%--------------------------------------------------------------------------
% Two-country Devereux-Sutherland (2011) model: Impulse responses.
% Version with real wealth, gross rates, and prices as states
% It replicates Figure 1, panel (b) of Appendix A.
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
gap = 0;

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d gap];
%make eta
eta=[zeros(4,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];
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
approx = 3; %order of approximation for dynamics

%color of plots?
my_color = 'b';
%my_color = 'r--';
%my_color = 'g';
%my_color = 'k';

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------

% First, solve for the SSS
sss_sol = compute_sss_ds_wealth_5p(model,params,M,eps_ind,[0 z0])
wh1 = sss_sol(1); %real holdings
z = sss_sol(2); %asset price

%Recalculate DSS of auxiliary model
gap = z/d0 - betta; %controls betta in aux. model
R = 1/(betta+gap); 
params(1) = wh1/z;
params(12) = gap;

%Deterministic steady-state
nxss=[wh1;0;z;z;zeros(6,1)]; %for states
nyss = [zeros(2,1);R;R;wh1/z;-wh1/z;z;z]; %for controls

mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1;
derivs1.hx(eps_ind,eps_ind) = 1; % evaluate at the model of interest (epsilon=1)

x1 = dr_ht(derivs1,nxss,approx,x0-nxss);
y1 = dr_gt(derivs1,nyss,approx,x0-nxss);
x1(eps_ind-1) = 0.01;

T0 = 1;
T = 11;
Tr = 10;
time_h = 2;
%make innovations
innovations = zeros(4,(T0 + (T-1))); % shocks from period 2 to T

%use Levintal's function simul.m to simulate the economy:
[yt,xt]=simul(x1,innovations,nyss,nxss,eta,derivs1,approx,0,model);

yt = [y1 yt];
xt = [x1 xt];

%make variables of interest
zht = yt(7,1:end-1);
zft = yt(8,1:end-1);
aht = yt(5,1:end-1);
aft = yt(6,1:end-1);
wealth = xt(2,2:end);

%Plot the results
time = 0:(Tr);

subplot(1,2,2);
plot(time,aht(T0:T0+Tr),my_color,'LineWidth',2);
title('(b) Portolio dynamics');

