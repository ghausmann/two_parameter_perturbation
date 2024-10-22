function y = compute_sss_ds(M_,options_,oo_,eps_ind,guess)
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% This function approximates the SSS of the Home asset
%
% Copyright (C) 2024 Guillermo Hausmann Guil

% Choose a second-order solution (enough for this model)
options_.order=2;
options_.k_order_solver=0;

%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_ds(x,M_, options_, oo_,eps_ind);

y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

