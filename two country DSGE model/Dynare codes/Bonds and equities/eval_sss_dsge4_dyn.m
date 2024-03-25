function y = eval_sss_dsge4_dyn(x, M_, options_, oo_, eps_ind)
% Two-country DSGE model with bonds and equities
% Baseline calibration with symmetric countries.
% This function implements the SSS algorithm, using a second-order
% approximation (enough for the SSS of this model).
%
% Inputs are x (vector of SSS values), Dynare's structs, and eps_ind (the
% index of the epsilon variable)
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of bonds and
%equities.
a1 = x(1);
b1 = x(2);

Sh1 = a1;
Sf1 = 1-a1;
Bh1 = b1;
Bf1 = -b1;

M_.params(1) = a1;
M_.params(2) = 1-a1;
M_.params(3) = b1;
M_.params(4) = -b1;

oo_.steady_state(3) = Bh1;
oo_.steady_state(4) = Bf1;
oo_.steady_state(5) = Sh1;
oo_.steady_state(6) = Sf1;

yss = oo_.steady_state;

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices of derivatives H1,..,Hk .
%The external algorithm is Dynare's MATLAB function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0 = yss(3:13,1);
x0(eps_ind) = 1; %This imposes epsilon=1

%This function evaluates the k-order decision rule h_k.
x1 = dr_yt(mdr,yss,2,x0-yss(3:13,1),zeros(6,1));
eq1 =  x1(3)-b1;
eq2 = x1(5)-a1;

y = [eq1;eq2];