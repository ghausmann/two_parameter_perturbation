%--------------------------------------------------------------------------
% Multi-asset DSGE model with bonds and equities: Comparative statics
%
% This script replicates Figure 7 in Appendix E, from the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;

%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
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
%Comparative statics: External wealth with asymmetric countries
%--------------------------------------------------------------------------

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
    % this function computes the SSS of bonds and equities when countries
    % are aysmmetric
    sol_gap(t,:) = compute_sss_dsge4_asym(M_,options_,oo_,eps_ind,[a (1-a) b -b]);
    
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
    % Appendix E.2).
    myfun = @(z)dsge_ss_y4(z,[phi alpa md betta],[a1 as1],[b1 bs1]);
    pf0 = fzero(myfun,1);
    
    P0 = ( alpa + (1-alpa)*pf0^(1-phi) )^(1/(1-phi));
    Ps0 = ( (1-alpa) + alpa*pf0^(1-phi) )^(1/(1-phi));
    C0 = (1/P0)*( (1-md) + a1*md + as1*pf0*md + (1-betta)*(b1*P0 + bs1*Ps0) );
    Cs0 = (1/Ps0)*( (1-md)*pf0 + (1-a1)*md + (1-as1)*pf0*md  -(1-betta)*(b1*P0 + bs1*Ps0) );
    
    %DSS asset prices
    zSf0 = pf0*(betta/(1-betta))*md;
    zBh0 = betta*P0;
    zBf0 = betta*Ps0;
    
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
    
    %SSS asset prices
    zBh_sss = zBh0*exp(y1(1));
    zBf_sss = zBf0*exp(y1(2));
    zSh_sss = zSh0*exp(y1(15));
    zSf_sss = zSf0*exp(y1(16));
    
    %SSS external positions
    ext_ap_gap(t) = (Sf0.*zSf_sss);
    ext_lp_gap(t) = ((1-Sh0).*zSh_sss);
    ext_ab_gap(t) = (Bh0.*zBh_sss);
    ext_lb_gap(t) = ((-Bf0).*zBf_sss);
    
end
toc

figure;
plot(gap_s,ext_ap_gap,'b');
hold on;
plot(gap_s,ext_ab_gap,'b--');
plot(gap_s,ext_lp_gap,'rx-');
plot(gap_s,ext_lb_gap,'rx--');
legend('equity assets','debt assets','equity liabilities','debt liabilities');
title('External wealth with asymmetric countries');
