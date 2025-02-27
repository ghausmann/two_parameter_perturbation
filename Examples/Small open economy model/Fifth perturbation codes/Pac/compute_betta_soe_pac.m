function y = compute_betta_soe_pac(P,b1,model,params,M,nxss0,nyss,eps_ind,approx,guess)
% This function calibrates the value of pssi2 (and indirectly the discount
% factor at the model of interest) consistent with a target SSS of bonds b_bar,
% with PAC as the auxiliary model.

%Evaluate the matrix eta:
ty = P(1);
tz = P(2);
params(6) = ty;
params(7) = tz;
eta=[zeros(1,2);ty 0;0 tz;zeros(1,2)];

%Call fsolve to find a candidate psi2 that zeroes the residual calibration
%function:
my_evalf= @(x)calib_sss_soe_pac(x,b1,model,params,M,eta,nxss0,nyss,eps_ind,approx);

y = fsolve(my_evalf,guess,optimset('TolFun',1e-15,'Display','off'));
