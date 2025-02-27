function y = compute_sss_soe_pac(P,model,params,M,nxss0,nyss,eps_ind,approx,guess)
% This function approximates the SSS of bonds b_bar of the SOE model with
% PAC as the auxiliary model.

%Evaluate the matrix eta:
ty = P(1);
tz = P(2);
params(6) = ty;
params(7) = tz;
eta=[zeros(1,2);ty 0;0 tz;zeros(1,2)];

%Call fsolve to find a candidate b_bar that zeroes the residual SSS
%function:
my_evalf= @(x)eval_sss_soe_pac(x,model,params,M,eta,nxss0,nyss,eps_ind,approx);

y = fsolve(my_evalf,guess,optimset('TolFun',1e-15,'Display','off'));

