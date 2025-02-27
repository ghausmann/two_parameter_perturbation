function y = compute_sss_growth(sz,model,params,M,eps_ind,approx,guess)
% This function approximates the SSS of capital k_bar in the neoclassical growth model.

%Evaluate the matrix eta:
eta=[0;sz;0];
params(6) = sz;

%Call fsolve to find a candidate b_bar that zeroes the residual SSS
%function:
fmy_evalf= @(x)eval_sss_growth(x,model,params,M,eta,eps_ind,approx);

y = fsolve(fmy_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
