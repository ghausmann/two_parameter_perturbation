function y = compute_calib_dsge4_dyn(a1,M_,options_,oo_,eps_ind,guess)
% This function solves for the std. of preference shocks and the SSS of
% bonds, such that the outcome is consistent with a SSS target for
% equities.
% Inputs are a1 a1 (target SSS of equities), Dynare's structs, eps_ind (index of the epsilon variable), and
% a initial guess.
% The oputput is the SSS solution.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

% Choose a second-order solution (enough for this model)
options_.order=2;
options_.k_order_solver=0;

% This is the residual function described in Section 2.4 of the paper, in
% hybrid calibration mode.
my_evalf= @(x)calib_sss_dsge4_dyn(x,a1,M_, options_, oo_,eps_ind);
% Use a non-linear solver to find the SSS.
mysol = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

y = mysol;