function y = compute_calib_dsge4_5p(a,model,params,M,eps_ind,guess)
% This function approximates the SSS of equities, for the
% baseline calibration with symmetric countries.

%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)calib_sss_dsge4_5p(x,a,model,params,M,eps_ind);
             
y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));

