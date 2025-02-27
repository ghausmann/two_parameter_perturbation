%--------------------------------------------------------------------------
% Multi-asset DSGE model with bonds and equities: Comparative statics
%
% This script replicates Figure 6 from the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;

% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
a = 0.83; %Just a initial value
af = 1-a;
b = 0; %Just a initial value
bf = -b;
betta = 0.96; % standard value
gama = 2; % standard value
pssi = 0.0001; % improves accuracy
pssi2 = 0.1; % improves accuracy
phi = 2; % Imbs and Mejean (2015)
alpa = 0.85; % U.S. import share
md = 0.036; %FRED data
d0 = log(md/(1-md));
rho_d = 0.42; %FRED data
td = 0.059; %FRED data
rho_y = 0.51; %FRED data
ty = 0.018; %FRED data
rho_q = 0.46; % within range of rho_d and rho_y
tq = 0.0067; % comparative statics
rho_eps = 1; % makes perturbation variable constant over time
kappa = 0.007; % match serial correlation of U.S. trade balance
%Deterministic steady-state values
C0 = 1;
Cs0 = 1;
P0 = 1;
Ps0 = 1;
pf0 = 1;
zSh0 = (betta/(1-betta))*md;
zSf0 = (betta/(1-betta))*md;
zBh0 = betta;
zBf0 = betta;
gap = 0; %symmetric case (gap>0 introduces assymetric risk aversion)
uy_uys_corr = 0.68; % Corsetti, Dedeloa and Leduc (2008)
uy_ud_corr = 0.12; %FRED data

params=[a af b bf betta gama pssi phi alpa md d0 rho_y ty rho_d td rho_q tq rho_eps kappa C0 Cs0 P0 Ps0 pf0 zSh0 zSf0 zBh0 zBf0 gap pssi2];
eta=[zeros(4,6);ty 0 0 0 0 0;0 ty 0 0 0 0;0 0 td 0 0 0;0 0 0 td 0 0;0 0 0 0 tq 0;0 0 0 0 0 tq;0 0 0 0 0 0];
%Moments
% Compute the cross moments:
n_e=6; % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 uy_uys_corr uy_ud_corr 0 0 0;uy_uys_corr 1 0 uy_ud_corr 0 0;uy_ud_corr 0 1 0 0 0;0 uy_ud_corr 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
M.M2 = Sigma;

eps_ind = 11; %index of epsilon in state vector
algo='gensylv';

%--------------------------------------------------------------------------
%Comparative statics: trade and financial linkages
%--------------------------------------------------------------------------

import_s = 0.1:0.01:0.45;
trade_open = 2*import_s;
l_is = length(import_s);
sol_imp = zeros(l_is,2);
ext_eq_imp = zeros(l_is,1);
ext_debt_imp = zeros(l_is,1);
share_eq_imp = zeros(l_is,1);
eq_home_bias = zeros(l_is,1);

for t=1:l_is
    
    alpa_t = 1-import_s(t);
    params(9) = alpa_t;
    % this function computes the SSS of bonds and equities
    sol_imp(t,:) = compute_sss_dsge4_5p(model,params,M,eps_ind,[1 1]);
                    
    a1 = sol_imp(t,1);
    b1 = sol_imp(t,2);
    
    % Recalculate DSS of auxiliary model
    Sh0 = a1;
    Sf0 = 1-a1;
    Bh0 = b1;
    Bf0 = -b1;
    
    params(1) = a1;
    params(2) = 1-a1;
    params(3) = b1;
    params(4) = -b1;
    
    nxss=[Bh0;Bf0;Sh0;Sf0;zeros(7,1)];
    nyss = zeros(11,1);
    
    %mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
    derivs1=solve_dsge(model,params,M,eta,nxss,nyss,2,algo);
    
    %Compute SSS of other variables
    x0=nxss; % start at the steady state
    x0(eps_ind) = 1; % evaluate at the model of interest (epsilon=1)
    ysss = dr_gt(derivs1,nyss,2,(x0-nxss));
    
    %SSS asset prices
    zBh_sss = zBh0*exp(ysss(3));
    zSf_sss = zSf0*exp(ysss(6));
    
    %SSS external positions
    ext_eq_imp(t) = (Sf0*zSf_sss);
    ext_debt_imp(t) = (Bh0*zBh_sss);
    
end

equity_open = 2*ext_eq_imp;
debt_open = 2*ext_debt_imp;
financial_open = equity_open + debt_open;


figure;
plot(trade_open,financial_open);
hold on;
plot(trade_open,equity_open,'g--');
plot(trade_open,debt_open,'r:');
legend('financial openness','equity openness','debt openness');
xlabel('Trade openness');