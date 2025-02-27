function y = compute_sss_dsge4(M_,options_,oo_,eps_ind,guess)
% This function approximates the SSS of bonds and equities, for the
% baseline calibration with symmetric countries.
% Inputs are Dynare's structs, eps_ind (index of the epsilon variable), and
% a initial guess.
% The oputput is the SSS.

% Choose a second-order solution (enough for this model)
options_.order=2;
options_.k_order_solver=0;

% This is the residual function described in Section 2.4 of the paper.
my_evalf= @(x)eval_sss_dsge4(x,M_, options_, oo_,eps_ind);
% Use a non-linear solver to find the SSS.
mysol = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

y = mysol;