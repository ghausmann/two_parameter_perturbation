function y = calib_sss_bond2(x,model,params,M,eta,eps_ind,approx,b_target)
% One bond, two-country model
% This function implements the SSS algorithm in calibration mode.
%
% The key input is x (the gap between countries' impatience). Other inputs
% are data file 'model', vector of parameters 'params', matrices of moments
% M, the 'eta' matrix, the index of epsilon in vector xt, the approximation
% order, and the target for NFA.

%STEP 1: Given b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.
gap = x;
b_bar = b_target;
betta = params(1);
C0h = 1 + (1-betta)*b_bar;
C0f = 1 - (1-betta)*b_bar;
R0 = 1/betta;

params(6) = b_bar;
params(7) = C0h;
params(8) = C0f;
params(9) = gap;

nxss = [b_bar;0;0;0];
nyss = [R0;1;1];

%STEP 2: Using the new DSS and other function inputs, call an external algorithm
%to compute the matrices H1,..,Hk of derivatives.
%The external algorithm is Levintal's function solve_dsge.m
algo='gensylv'; %algorithm to compute first-order derivatives
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0=nxss;
if eps_ind>0
    x0(eps_ind) = 1; % impose the model of interest
    derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon follows a unit-root process.
end
%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
x1=dr_ht(derivs1,nxss,approx,(x0-nxss));

%Return h1(0,1,1)-b_bar
y = x1(1)-b_bar;