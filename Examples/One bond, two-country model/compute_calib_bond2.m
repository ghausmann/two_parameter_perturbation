function y = compute_calib_bond2(model,params,M,eta,eps_ind,approx,guess,b_target)
% This function uses the SSS algorithm in calibration mode to compute the
% gap between countries' impatience consistent with a SSS NFA target for the
% Home country.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

my_evalf= @(x)calib_sss_bond2(x,model,params,M,eta,eps_ind,approx,b_target);

y = fsolve(my_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
