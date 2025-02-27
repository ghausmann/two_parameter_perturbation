function y = eval_sss_soe_uzawa(b_bar,model,params,M,eta,nxss0,nyss,eps_ind,approx)
% Small open economy model
% Auxiliary model is Uzawa
% This function implements the SSS algorithm.
%
% The key input is b_bar.
% Other inputs are data file 'model', vector of parameters 'params', 
% matrices of moments M, the 'eta' matrix, the DSS of exogenous states
% 'nxss0', the DSS of controls 'nyss', the index of epsilon in vector xt,
% and the approximation order.

%STEP 1: Given b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.
R0 = params(8);
C01 = 1 + (1-(1/R0))*b_bar;
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
x0=nxss; 
if eps_ind>0
    x0(eps_ind) = 1; % impose the model of interest
    derivs1.hx(eps_ind,eps_ind) = 1;
end
%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
xt=dr_ht(derivs1,nxss,approx,(x0-nxss));

%Return h1(0,1,1)-b_bar 
y = xt(1)-b_bar;
