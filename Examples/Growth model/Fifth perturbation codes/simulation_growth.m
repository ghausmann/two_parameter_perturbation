% -------------------------------------------------------------------------
% Neoclassical growth model: Simulations
% -------------------------------------------------------------------------

clear
% Seed the pseudo-random number generator, uncomment for replication.
rng(0)
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose the calibration. Zero for standard, and one for extreme.
% Note: default is extreme calibration.
calibration = 1;
% load the global solution
load('my_k_poli');

% Choose parameter values:
betta = 0.99;
alpa = 0.35;
delta = 0.025;
rho_z = 0.95;
rho_eps = 0.999999; % very close (or equal) to 1
psi = 0; % just to provide a initial value

if calibration==0
    gama = 2;
    sz = 0.01;
elseif calibration==1
    gama = 20;
    sz = 0.025;
end

% deterministic steady-state
Kss_0 = (alpa/(1/betta - (1-delta)))^(1/(1-alpa));
% Make eta
eta=[0;sz;0];
% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.

% Choose algorithm for the Sylvester equation
algo='gensylv'; % Kamenik algorithm
approx = 3; % order of approximation
eps_ind = 3; % index of the epsilon variable

% -------------------------------------------------------------------------
% Standard perturbation
% -------------------------------------------------------------------------

% deterministic steady-state of the model (in logs)
k_ss1 = log(Kss_0); % capital
ey_ss1 = log(exp(k_ss1)^alpa); %output
ei_ss1 = log(delta*exp(k_ss1)); %investment
c_ss1 = log(exp(ey_ss1) - exp(ei_ss1)); %consumption
z_ss = 0; % log of TPF
eps_ss = 0; % epsilon perturbation variable

% collect deterministic steady-state
nxss1=[k_ss1;z_ss;eps_ss]; % states
nyss1=[c_ss1;ei_ss1;ey_ss1]; % controls

% Make a vector of parameter values:
params1=[betta gama alpa delta rho_z sz rho_eps psi];
% check deterministic steady-state
mycheck1 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss1(:);nyss1(:);nxss1(:);nxss1(:);params1(:)]));
% solve perturbation model
derivs1 = solve_dsge(model,params1,M,eta,nxss1,nyss1,approx,algo);

x0_1 = nxss1; % initial state for the simulation

% -------------------------------------------------------------------------
% Two-parameter perturbation
% -------------------------------------------------------------------------

% compute stochastic steady-state of the model of interest, for a given
% level of risk (auxiliary parameter psi will change within the function)
k2 = compute_sss_growth(sz,model,params1,M,eps_ind,approx,Kss_0);

% implied calibration of the auxiliary parameter:
psi = betta*(alpa/(k2^(1-alpa)) + (1-delta)) - 1;
% Make a new vector of parameter values:
params2 = params1;
params2(8) = psi; % update psi

% new deterministic steady-state of the auxiliary model (in logs)
k_ss2 = log(k2); % capital
ey_ss2 = log(exp(k_ss2)^alpa); %output
ei_ss2 = log(delta*exp(k_ss2)); %investment
c_ss2 = log(exp(ey_ss2) - exp(ei_ss2)); %consumption

% collect deterministic steady-state of the auxiliary model
nxss2=[k_ss2;z_ss;eps_ss]; % states
nyss2=[c_ss2;ei_ss2;ey_ss2]; % controls

% check deterministic steady-state of the auxiliary model
mycheck2 = double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss2(:);nyss2(:);nxss2(:);nxss2(:);params2(:)]));
% solve the two-parameter perturbation model
derivs2=solve_dsge(model,params2,M,eta,nxss2,nyss2,approx,algo);

% initial state for the simulation, same as in standard perturbation..
x0_2 = nxss1;
%..but evaluating at the model of interest..
x0_2(eps_ind) = 1;
derivs2.hx(eps_ind,eps_ind) = 1; % ..and imposing full unit-root behavior

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

% paths of capital for different solution methods
kt_global = zeros(T0+T,1); % global
kt_global(1) = exp(x0_2(1));

kt_local1 = zeros(T0+T,1); % standard perturbation
kt_local1(1) = exp(x0_1(1));

kt_local2 = zeros(T0+T,1); % two-parameter perturbation
kt_local2(1) = exp(x0_2(1));

% States for local solutions
mxt1_0 = x0_1; % standard perturbation
mxt2_0 = x0_2; % two-parameter perturbation

% implement simulation
for t=2:T0+T

    % update capital with global solution:
    kt_global(t) = interp1(k_grid,k_poli(my_simul(t),:),kt_global(t-1),'spline','extrap');

    % update exogenous state z:
    mxt1_0(2) = zz(my_simul(t));
    mxt2_0(2) = zz(my_simul(t));

    % update capital with local solutions:
    % standard perturbation
    mxt1_1 = dr_ht(derivs1,nxss1,approx,(mxt1_0-nxss1));
    kt_local1(t) = exp(mxt1_1(1));
    % two-parameter perturbation
    mxt2_1 = dr_ht(derivs2,nxss2,approx,(mxt2_0-nxss2));
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
