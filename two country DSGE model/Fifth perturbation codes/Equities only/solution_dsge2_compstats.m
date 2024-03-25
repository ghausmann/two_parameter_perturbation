%--------------------------------------------------------------------------
% Two-country DSGE model with bonds and equities: Comparative statics
% This script replicates (part of) Figure 1 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;
% Add folder 'files' to the search path:
addpath('files'); 
% Load the model:
load('model')

% Choose parameter values:
a = 0.5; %Just a initial value
betta = 0.96; %standard value
gama = 2; %standard value
pssi = 0.001; % improves accuracy
phi = 2; %Imbs and Mejean (2015)
alpa = 0.85; % U.S. import share
md = 0.036; % DSS dividend share (FRED data)
d0 = log(md/(1-md));
rho_d = 0.42; %FRED data
td = 0.059; %FRED data
rho_y = 0.51; %FRED data
ty = 0.018; %FRED data
rho_q = 0.46; % within range of rho_d and rho_y
tq = 0.0001; % comparative statics
rho_eps = 1; % makes perturbation variable constant over time
kappa = 0.007; % match serial correlation of U.S. trade balance
uy_uys_corr = 0.68; % Corsetti, Dedeloa and Leduc (2008)
uy_ud_corr = 0.12; %FRED data
%Deterministic steady-state values
C0 = 1;
Cs0 = 1;
P0 = 1;
Ps0 = 1;
pf0 = 1;
zSh0 = (betta/(1-betta))*md;
zSf0 = (betta/(1-betta))*md;

params=[a betta gama pssi phi alpa md d0 rho_y ty rho_d td rho_q tq rho_eps kappa C0 Cs0 P0 Ps0 pf0 zSh0 zSf0];

%Moments
% Compute the cross moments:
n_e=6; % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 uy_uys_corr uy_ud_corr 0 0 0;uy_uys_corr 1 0 uy_ud_corr 0 0;uy_ud_corr 0 1 0 0 0;0 uy_ud_corr 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]; 
M.M2 = Sigma;

eps_ind = 9; %index of epsilon in state vector

%--------------------------------------------------------------------------
%Comparative statics: std. redistributive shocks
%--------------------------------------------------------------------------

mtd = 0.025:0.005:0.1;
l_mtd = length(mtd);
sol_mtd = zeros(l_mtd,1);

for t=1:l_mtd
    
  params(12)=mtd(t); % update std. 
  % this function computes the SSS of equities
  sol_mtd(t) = compute_sss_dsge2(model,params,M,eps_ind,0.5);
      
end

%--------------------------------------------------------------------------
%Comparative statics: std. preference shocks
%--------------------------------------------------------------------------
params(12)=td; % set td back to its calibrated value

mtq = 0.0001:0.0005:0.01;
l_mtq = length(mtq);
sol_mtq = zeros(l_mtq,1);

for t=1:l_mtq
    
  params(14)=mtq(t); % update std.
  % this function computes the SSS of equities
  sol_mtq(t) = compute_sss_dsge2(model,params,M,eps_ind,0.5);
       
end

% Plot the results
figure;
subplot(1,2,1);
plot(mtd,sol_mtd(:,1),'r');
title('(a) As a function of redistributive risk');
xlabel('\eta_{d}');
subplot(1,2,2);
plot(mtq,sol_mtq(:,1),'r');
title('(b) As a function of preference risk');
xlabel('\eta_{q}');


