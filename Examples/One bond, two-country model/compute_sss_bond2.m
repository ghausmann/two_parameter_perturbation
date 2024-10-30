function y = compute_sss_bond2(model,params,M,eta,eps_ind,approx,guess)
% This function uses the SSS algorithm to compute the SSS of Home NFA in
% the one-bond, two-country model.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%Call fsolve to find a candidate b_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_bond2(x,model,params,M,eta,eps_ind,approx);

y = fsolve(my_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
