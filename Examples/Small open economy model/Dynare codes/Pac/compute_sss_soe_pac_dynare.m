function y = compute_sss_soe_pac_dynare(M_,options_,oo_,eps_ind,guess)
%This function approximates the SSS of bonds b_bar of the SOE model with
%PAC as the auxiliary model.
%Output is the SSS of bonds.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

options_.order=2; %Choose a second-order solution
options_.k_order_solver=0;

%Call fsolve to find a candidate DSS for bonds that zeroes the residual
%calibration function:
my_evalf= @(x)eval_sss_soe_pac_dynare(x,M_,options_,oo_,eps_ind);

y = fsolve(my_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
