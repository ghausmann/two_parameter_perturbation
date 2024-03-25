%--------------------------------------------------------------------------
% Two-country DSGE model with bonds and equities: Comparative statics
% This script replicates Figure 2 in the paper.
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

sq = 0.0067 -3.4135e-05; %Calibrated std. of preference shocks
M_.params(17) = sq;

eps_ind = 5; %index of perturbation variable

%--------------------------------------------------------------------------
%Comparative statics: import share
%--------------------------------------------------------------------------

import_s = 0.1:0.01:0.45;
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

%--------------------------------------------------------------------------
%Comparative statics: heterogeneous countries
%--------------------------------------------------------------------------
M_.params(9) = alpa; % set alpa back to its calibrated value

gap_s = 0:0.05:0.30;
l_g = length(gap_s);
sol_gap = zeros(l_g,4);
ext_ap_gap = zeros(l_g,1);
ext_lp_gap = zeros(l_g,1);
ext_ab_gap = zeros(l_g,1);
ext_lb_gap = zeros(l_g,1);

guess = [a (1-a) b -b];

tic
for t=1:l_g
    
    gap_t = -gap_s(t);
    M_.params(31) = gap_t;
    sol_gap(t,:) = compute_sss_dsge4_asym_dyn(M_,options_,oo_,eps_ind,[a (1-a) b -b]);
    
    guess = sol_gap(t,:);
    a1 =  sol_gap(t,1);
    as1 =  sol_gap(t,2);
    b1 =  sol_gap(t,3);
    bs1 =  sol_gap(t,4);
    
    Sh0 = a1;
    Sf0 = as1;
    Bh0 = b1;
    Bf0 = bs1;
    
    M_.params(1) = a1;
    M_.params(2) = as1;
    M_.params(3) = b1;
    M_.params(4) = bs1;
    
    % Solve for the DSS of the terms of trade with a non-linear solver (see
    % Appendix C).
    myfun = @(z)dsge_ss_y4(z,[phi alpa md betta],[a1 as1],[b1 bs1]);
    pf0 = fzero(myfun,1);
    
    P0 = ( alpa + (1-alpa)*pf0^(1-phi) )^(1/(1-phi));
    Ps1 = ( (1-alpa) + alpa*pf0^(1-phi) )^(1/(1-phi));
    C0 = (1/P0)*( (1-md) + a1*md + as1*pf0*md + (1-betta)*(b1*P0 + bs1*Ps1) );
    Cs0 = (1/Ps1)*( (1-md)*pf0 + (1-a1)*md + (1-as1)*pf0*md  -(1-betta)*(b1*P0 + bs1*Ps1) );
    
    %DSS asset prices
    zSf0 = pf0*(betta/(1-betta))*md;
    zBh0 = betta*P0;
    zBf0 = betta*Ps1;
    
    M_.params(20) = C0;
    M_.params(21) = Cs0;
    M_.params(22) = P0;
    M_.params(23) = Ps0;
    M_.params(24) = pf0;
    M_.params(26) = zSf0;
    M_.params(27) = zBh0;
    M_.params(28) = zBf0;
    
    oo_.steady_state(3) = Bh0;
    oo_.steady_state(4) = Bf0;
    oo_.steady_state(5) = Sh0;
    oo_.steady_state(6) = Sf0;
    
    %Compute the perturbation solution:
    %The external algorithm is Dynare matlab function resol.m
    [mdr, ~, ~, ~] = resol(0, M_, options_, oo_);
    yss = mdr.ys;
    
    %Compute SSS of other variables
    x0 = yss(3:13,1);
    x0(eps_ind) = 1;
    y1 = dr_yt(mdr,yss,2,x0-yss(3:13,1),zeros(6,1));
    
    zBh_sss = zBh0*exp(y1(1));
    zBf_sss = zBf0*exp(y1(2));
    zSh_sss = zSh0*exp(y1(15));
    zSf_sss = zSf0*exp(y1(16));
    
    ext_ap_gap(t) = (Sf0.*zSf_sss);
    ext_lp_gap(t) = ((1-Sh0).*zSh_sss);
    ext_ab_gap(t) = (Bh0.*zBh_sss);
    ext_lb_gap(t) = ((-Bf0).*zBf_sss);
    
end
toc

figure;
subplot(1,2,1);
plot(import_s,ext_eq_imp);
hold on;
plot(import_s,ext_debt_imp,'g--');
legend('equity','bonds');
title('(a) As a function of trade openness');
subplot(1,2,2);
plot(gap_s,ext_ap_gap,'b');
hold on;
plot(gap_s,ext_ab_gap,'b--');
plot(gap_s,ext_lp_gap,'rx');
plot(gap_s,ext_lb_gap,'rx--');
legend('equity assets','debt assets','equity liabilities','debt liabilities');
title('(b) Heterogeneous risk aversion');
