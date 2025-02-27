function y = compute_sss_growth_dyn(M_,options_,oo_,eps_ind,guess)
%This function approximates the SSS of capital k_bar in the neoclassical growth model.

%Call fsolve to find a candidate b_bar that zeroes the residual SSS
%function:
fmy_evalf= @(x)eval_sss_growth_dyn(x,M_, options_, oo_,eps_ind);

y = fsolve(fmy_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
