% -------------------------------------------------------------------------
% Neoclassical growth model: Comparative statics
% This script replicates the results in Figure 2 of Appendix D.2
% -------------------------------------------------------------------------
clear
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
betta = 0.99;
gama = 2;
alpa = 0.35;
delta = 0.025;
rho_z = 0.95;
sz = 0.01; % just to provide a initial value
rho_eps = 1; % very close (or equal) to 1
psi = 0; % just to provide a initial value

% deterministic steady-state
Kss_0 = (alpa/(1/betta - (1-delta)))^(1/(1-alpa));
% Make eta
eta=eval(eta);
% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.

approx = 3;
eps_ind = 3; %index of the epsilon variable

msz = 0.0025:0.0025:0.03; %grid for std. of innovations
lsz = length(msz);

% -------------------------------------------------------------------------
% baseline calibration: gamma = 2
% -------------------------------------------------------------------------
gama1 = 2;
params1=[betta gama1 alpa delta rho_z sz rho_eps psi];

% set a loop to compute stochastic steady-states for different levels of
% risk (as measured by sz).

my_ksss1_p1 = zeros(lsz,1);
my_ksss1_p2 = zeros(lsz,1);

for t=1:lsz

    % computes the stochastic steady-state with standard perturbation
    % by iterating the states until convergence
    my_ksss1_p1(t) = iter_sss_growth(msz(t),model,params1,M,approx);
    % computes the stochastic steady-state with two-parameter perturbation
    % by solving a non-linear equation
    my_ksss1_p2(t) = compute_sss_growth(msz(t),model,params1,M,eps_ind,approx,Kss_0);

end

% -------------------------------------------------------------------------
% extreme calibration: gamma = 20
% -------------------------------------------------------------------------
gama2 = 20;
params2=[betta gama2 alpa delta rho_z sz rho_eps psi];

my_ksss2_p1 = zeros(lsz,1);
my_ksss2_p2 = zeros(lsz,1);

for t=1:lsz

    my_ksss2_p1(t) = iter_sss_growth(msz(t),model,params2,M,approx);
    my_ksss2_p2(t) = compute_sss_growth(msz(t),model,params2,M,eps_ind,approx,Kss_0);

end

% load the comparative statics from the global solution
load('my_comp_sss_global');

% plot the results
figure;
subplot(1,2,1);
plot(msz,my_ksss1_g);
hold on;
plot(msz,my_ksss1_p1,'g--');
plot(msz,my_ksss1_p2,'r-.');

subplot(1,2,2);
plot(msz,my_ksss2_g);
hold on;
plot(msz,my_ksss2_p1,'g--');
plot(msz,my_ksss2_p2,'r-.');
