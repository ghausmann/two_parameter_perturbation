function y = compute_sss_ds_wealth_5p(model,params,M,eps_ind,guess)
% Devereux-Sutherland (2011) model
% Version with real wealth, gross rates, and prices as states
% This function approximates the SSS of equities, for the
% baseline calibration with symmetric countries.


%Evaluate the matrix eta:
tk = params(6);
tl = params(7);
% Make eta
eta=[zeros(4,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];
%Call fsolve to find a candidate a_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_ds_wealth_5p(x,model,params,M,eta,eps_ind);

y = fsolve(my_evalf,guess,optimset('Display','off','TolFun',1e-12));
