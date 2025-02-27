
function y = compute_sss_dsge2(M_,options_,oo_,eps_ind,guess)
% This function approximates the SSS of equities, for the
% baseline calibration with symmetric countries.

% Choose a second-order solution (enough for this model)
options_.order=2;
options_.k_order_solver=0;

%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_dsge2(x,M_, options_, oo_,eps_ind);

y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

