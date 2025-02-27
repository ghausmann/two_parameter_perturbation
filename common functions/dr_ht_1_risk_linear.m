function y = dr_ht_1_risk_linear(jac_hx_rl,hss_rl,x)
%This function evaluates linear decision rules for state variables of a DSGE
%model solved with the algorithm by Levintal (2017).
%
%Inputs are:
%jac_gx_rl: nx-by-nx matrix of coefficients.
%hss_rl: nx-by-1 vector of zero-order terms.
%x: nx-by-N vector of state variables, where N>=1, and each column is a particular
%realization of states.
%
%The output is a nx-by-N matrix, where each column represents time t+1
%state variables (excluding future innovations), under a particular
%realization of today's states.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[~,s2] = size(x);

x = [x;ones(1,s2)];
hx = jac_hx_rl;


y = hss_rl +  hx*x;

