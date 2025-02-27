function g = dr_gt_1_risk_linear(jac_gx_rl,gss_rl,x)
%This function evaluates linear  decision rules for control variables of a DSGE
%model solved with the algorithm by Levintal (2017).
%
%Inputs are:
%jac_gx_rl: ny-by-nx matrix of coefficients.
%gss_rl: ny-by-1 vector of zero-order terms.
%x: nx-by-N vector of state variables, where N>=1, and each column is a particular
%realization of states.
%
%The output is a ny-by-N matrix, where each column represents
%time t control variables under a particular realization of states.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[~,s2] = size(x);

x = [x;ones(1,s2)];

gx = jac_gx_rl;

g = gss_rl +  gx*x;


