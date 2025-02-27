% -------------------------------------------------------------------------
% Neoclassical growth model: Simulations
% This script replicates the results in Table 1 of Appendix D.2
% -------------------------------------------------------------------------

%NOTE: make sure to add the master folder "Growth model" to the search
%path, including subfolders.

clear
% Seed the pseudo-random number generator, uncomment for replication.
rng(0)
%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
%Load the pre-processing data
load('my_growth_model.mat');

% Choose the calibration. Zero for standard, and one for extreme.
% Default is extreme calibration.
calibration = 1;
% load the global solution
load('my_k_poli');

if calibration==0
    gama = 2;
    sz = 0.01;
elseif calibration==1
    gama = 20;
    sz = 0.025;
end

M_.params(2) = gama;
M_.params(6) = sz;

% index of the epsilon variable 
eps_ind_y = 4; %in the vector yt (Dynare's endogenous variables)
eps_ind_x = 3; %in the vector xt (state variables)

% deterministic steady-state
Kss_0 = (alpa/(1/betta - (1-delta)))^(1/(1-alpa));

% -------------------------------------------------------------------------
% Standard perturbation
% -------------------------------------------------------------------------

% Solve perturbation model, using Dynare's function "resol.m"
[mdr1, ~, M_1, ~] = resol(0, M_, options_, oo_);

y0_1 = mdr1.ys; % initial state for the simulation

% -------------------------------------------------------------------------
% Two-parameter perturbation
% -------------------------------------------------------------------------

% Compute the stochastic steady-state of the model of interest (the auxiliary
% parameter psi will change within the function).
% will take a few seconds because it uses Dynare's third-order solver..
tic
k2 = compute_sss_growth_dyn(M_,options_,oo_,eps_ind_x,1.4*Kss_0);
toc
% implied calibration of the auxiliary parameter:
psi = betta*(alpa/(k2^(1-alpa)) + (1-delta)) - 1;
% new deterministic steady-state of the auxiliary model (in logs)
kss2 = log(k2); % capital
eyss2 = log(exp(kss2)^alpa); %output
eiss2 = log(delta*exp(kss2)); %investment
css2 = log(exp(eyss2) - exp(eiss2)); %consumption

% collect steady-state values
yss2 = [eiss2;kss2;0;0;css2;eyss2];
oo_.steady_state = yss2;
%change auxiliary parameter value to the implied one.
M_.params(8) = psi; 

% solve perturbation model
[mdr2, ~, M_2, ~] = resol(0, M_, options_, oo_);

% initial state for the simulation, same as in standard perturbation..
y0_2 = y0_1;
%..but evaluating at the model of interest
y0_2(eps_ind_y) = 1;

% -------------------------------------------------------------------------
% Stochastic simulation
% -------------------------------------------------------------------------

T0 = 1000; % burn-in periods
T = 10000; % number of periods

% Generate a stochastic simulation for z, using the global approach (a
% markov chain over 11 grid-points for z generated with the Rouwenhorst
% method)
my_rand = rand([T0+T,1]);
my_simul = hitm_s(Pr',my_rand);

% path of capital for different solution methods
kt_global = zeros(T0+T,1); % global
kt_global(1) = exp(y0_1(2));

kt_local1 = zeros(T0+T,1); % standard perturbation
kt_local1(1) = exp(y0_1(2));

kt_local2 = zeros(T0+T,1); % two-parameter perturbation
kt_local2(1) = exp(y0_1(2));

% state variables for local solutions
mxt1_0 = y0_1(2:4,1); % standard perturbation
mxt2_0 = y0_2(2:4,1); % two-parameter perturbation

% implement simulation
for t=2:T0+T

    % update capital with global solution:
    kt_global(t) = interp1(k_grid,k_poli(my_simul(t),:),kt_global(t-1),'spline','extrap');

    % update exogenous state z
    % The trick: impose pre-determined z(t-1)=0, and just update the current shock, scaled
    % by 1/sz. This is equivalent to make z(t) the current exogenous state. 
    mxt1_0(2) = 0;
    mxt2_0(2) = 0;
    ut = zz(my_simul(t))/sz;
    
    % update capital with local solutions, using my own function "dr_yt.m"
    % standard perturbation
    myt1_1 = dr_yt(mdr1,mdr1.ys,3,mxt1_0-mdr1.ys(2:4,1),ut);
    mxt1_1 = myt1_1(2:4,1); %pick the states only
    kt_local1(t) = exp(mxt1_1(1));
    % two-parameter perturbation
    myt2_1 = dr_yt(mdr2,mdr2.ys,3,mxt2_0-mdr2.ys(2:4,1),ut);
    mxt2_1 = myt2_1(2:4,1); %pick the states only
    kt_local2(t) = exp(mxt2_1(1));

    mxt1_0 = mxt1_1;
    mxt2_0 = mxt2_1;

end

%Plot the time-series of capital under different solution methods
figure;plot(kt_global(T0+1:T0+2000));
hold on;
plot(kt_local1(T0+1:T0+2000),'g');
plot(kt_local2(T0+1:T0+2000),'r');
title('time series of capital');

% compute percent errors (relative to global solution) and report stats
percent_error1 = 100*abs(kt_local1(T0+1:end)-kt_global(T0+1:end))./kt_global(T0+1:end);
stats1 = [mean(percent_error1) max(percent_error1)]

percent_error2 = 100*abs(kt_local2(T0+1:end)-kt_global(T0+1:end))./kt_global(T0+1:end);
stats2 = [mean(percent_error2) max(percent_error2)]

%plot percent errors
figure;plot(percent_error1,'g');hold on;plot(percent_error2,'r');
title('percent errors');
