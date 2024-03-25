function y = compute_sss_dsge4(model,params,M,eps_ind,guess)
% This function approximates the SSS of equities, for the
% baseline calibration with symmetric countries.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%Evaluate the matrix eta:
ty = params(13);
td = params(15);
tq = params(17);
% Make eta
eta=[zeros(4,6);ty 0 0 0 0 0;0 ty 0 0 0 0;0 0 td 0 0 0;0 0 0 td 0 0;0 0 0 0 tq 0;0 0 0 0 0 tq;0 0 0 0 0 0];

%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_dsge4(x,model,params,M,eta,eps_ind);

y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

