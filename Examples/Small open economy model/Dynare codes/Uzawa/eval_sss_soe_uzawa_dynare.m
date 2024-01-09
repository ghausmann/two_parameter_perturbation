function y = eval_sss_soe_uzawa_dynare(b_bar,M_, options_, oo_,eps_ind)
% Small open economy model
% Auxiliary model is Uzawa.
% This function implements the algorithm described in Section 2.4 of the
% paper.
%Inputs are b_bar (SSS of bonds), structs generated by Dynare (M_,
%options_, and oo_), and the index of the epsilon variable in the vector yt
%of endogenous variables.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.
R0 = M_.params(8);
C01 = 1 + (1-(1/R0))*b_bar;
M_.params(9) = b_bar;
M_.params(10) = C01;
yss = [b_bar;0;0;0;1];
oo_.steady_state = yss;

%STEP 2: Using the new DSS and other function inputs, call an external algorithm
%to compute the matrices H1,..,Hk of derivatives.
%The external algorithm is Dynare matlab function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest epsilon=sigma=1.
x0 = yss(1:4,1);
x0(eps_ind) = 1;

%Compute the Taylor-series policy rule for bonds
%uses my the function dr_yt.m
x1 = dr_yt(mdr,yss,2,x0-yss(1:4,1),[0;0]);

y = x1(1) - b_bar;
