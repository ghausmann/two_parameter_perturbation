function y = eval_sss_dsge2_dyn(a1,M_, options_, oo_,eps_ind)
% Two-country DSGE model with equities only
% Baseline calibration with symmetric countries.
% This function implements the algorithm described in Section 2.4 of the
% paper.
% Inputs are x (vector of SSS values), Dynare's structs, and eps_ind (the
% index of the epsilon variable)
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of equities.
Sh1 = a1;
Sf1 = 1-a1;

M_.params(1) = a1;

oo_.steady_state(1) = Sh1;
oo_.steady_state(2) = Sf1;

yss = oo_.steady_state;

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices of derivatives H1,..,Hk .
%The external algorithm is Dynare's MATLAB function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0 = yss(1:9,1);
x0(eps_ind) = 1; %This imposes epsilon=1

%This function evaluates the k-order decision rule h_k.
x1 = dr_yt(mdr,yss,2,x0-yss(1:9,1),zeros(6,1));

y = x1(1) - a1;