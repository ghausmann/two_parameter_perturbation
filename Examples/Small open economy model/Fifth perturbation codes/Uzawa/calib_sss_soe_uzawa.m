function y = calib_sss_soe_uzawa(psi2,b_bar,model,params,M,eta,nxss0,nyss,eps_ind,approx)
% Small open economy model
% Auxiliary model is Uzawa
%
% This function implements the SSS algorithm in calibration mode.
% We pin down pssi2 (controls the discount factor of the model of interest)
% consistent with a fixed target for SSS bonds b_bar.
% Other inputs are data file 'model', the vector of parameters 'params', the
% matrices of moments M, the 'eta' matrix, the DSS of exogenous states
% 'nxss0', the DSS of controls 'nyss', the index of epsilon in vector xt,
% and the approximation order.


%STEP 1: Given psi2 and b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.

%Assign new value to psi2
params(12) = psi2;
%Recalculate DSS
R0 = params(8);
b_bar1 = b_bar;
C01 = 1 + (1-(1/R0))*b_bar1;
params(9) = b_bar;
params(10) = C01;
nxss=[b_bar;nxss0];

%STEP 2: Using the new DSS and other function inputs, call an external algorithm
%to compute the matrices H1,..,Hk of derivatives.
%The external algorithm is Levintal's function solve_dsge.m
algo='gensylv'; %algorithm to compute first-order derivatives
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest epsilon=sigma=1.
x0=nxss; % start at the steady state
if eps_ind>0
    x0(eps_ind) = 1; % impose the model of interest
    derivs1.hx(eps_ind,eps_ind) = 1;
end
%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
xt=dr_ht(derivs1,0,approx,(x0-nxss));

%Return h1(0,1,1)-b_bar
y = xt(1);