function y = eval_sss_growth_dyn(k_bar,M_, options_, oo_,eps_ind)
%Neoclassical growth model
%This function implements the SSS algorithm of the paper. The key input is
%k_bar (DSS target for capital) Other inputs are structs generated by
%Dynare (M_,options_, and oo_), and the index of the epsilon variable in
%the vector yt of endogenous variables.
%The order of approximation is fixed at k=3.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given k_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.

betta = M_.params(1);
alpa = M_.params(3);
delta = M_.params(4);

k1 = k_bar(1); % approximation point for capital, in levels
psi = betta*(alpa/(k1^(1-alpa)) + (1-delta)) - 1; %implied value of psi, given k1

% implied steady-state (in logs)
kss = log(k1);
eyss = log(exp(kss)^alpa);
eiss = log(delta*exp(kss));
css = log(exp(eyss) - exp(eiss));
zss = 0;
epsss = 0;

% collect steady-state values
yss = [eiss;kss;zss;epsss;css;eyss];
oo_.steady_state = yss;

M_.params(8) = psi; %change auxiliary parameter value to the implied one.

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices H1,..,Hk of derivatives.
%The external algorithm is Dynare's MATLAB function resol.m
% options_.order=2;
% options_.k_order_solver=0;
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest epsilon=sigma=1.
x0 = yss(2:4,1); % start at the steady state
x0(eps_ind) = 1; % impose the model of interest (eps=1)

%Compute the Taylor-series policy rule for bonds
%uses the function dr_yt.m
y1 = dr_yt(mdr,yss,3,x0-yss(2:4,1),0);
x1 = y1(2:4,1); %pick the states only

%Return h1(0,1,1)- log(k_bar)
y = x1(1) - log(k1);