%--------------------------------------------------------------------------
% Two-country DSGE model with bonds and equities: Comparative statics
% This script replicates Figure 5 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;

%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.2\matlab');
%Load the pre-processing data
load('my_dsge4_model.mat');

% Choose a second-order solution (enough for this model)
options_sss = options_;
options_sss.order=2;
options_sss.k_order_solver=0;

tq = 0.0067; %Calibrated std. of preference shocks
M_.params(17) = tq;

eps_ind = 5; %index of perturbation variable

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

tic
for t=1:l_is
    
    alpa_t = 1-import_s(t);
    M_.params(9) = alpa_t;
    % this function computes the SSS of bonds and equities
    sol_imp(t,:) = compute_sss_dsge4_dyn(M_,options_sss,oo_,eps_ind,[1 1]);
    
    a1 = sol_imp(t,1);
    b1 = sol_imp(t,2);
    
    % Recalculate DSS of auxiliary model
    Sh0 = a1;
    Sf0 = 1-a1;
    Bh0 = b1;
    Bf0 = -b1;
    
    M_.params(1) = a1;
    M_.params(2) = 1-a1;
    M_.params(3) = b1;
    M_.params(4) = -b1;
    
    oo_.steady_state(3) = Bh0;
    oo_.steady_state(4) = Bf0;
    oo_.steady_state(5) = Sh0;
    oo_.steady_state(6) = Sf0;
    yss = oo_.steady_state;
    
    %Compute the perturbation solution:
    %The external algorithm is Dynare matlab function resol.m
    [mdr, ~, ~, ~] = resol(0, M_, options_sss, oo_);
    
    %Compute SSS of other variables
    x0 = yss(3:13,1);
    x0(eps_ind) = 1; %evaluate at the model of interest
    y1 = dr_yt(mdr,yss,2,x0-yss(3:13,1),zeros(6,1));
    
    %SSS asset prices
    zBh_sss = zBh0*exp(y1(1));
    zSh_sss = zSh0*exp(y1(15));
    zSf_sss = zSf0*exp(y1(16));
    
    %SSS external positions
    ext_eq_imp(t) = (Sf0*zSf_sss);
    ext_debt_imp(t) = (Bh0*zBh_sss);
    
end
toc

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
