function y = calib_sss_dsge4(x,a1,model,params,M,eps_ind)
% Two-country DSGE model with bonds and equities
% This function implements the SSS algorithm in hybrid calibration mode.
%
% Inputs are the vector x (first component is the std. of preference
% shocks, and second is the SSS of bonds), a1 (SSS target for equities),
%data file 'model', vector of parameters 'params', matrices of moments M,
%and eps_ind (the index of the epsilon variable).
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we need to recalculate the DSS of bonds and
%equities, and update the eta matrix.
mtq = x(1);
b1 = x(2);

params(17) = mtq;
ty = params(13);
td = params(15);

Sh1 = a1;
Sf1 = 1-a1;
Bh1 = b1;
Bf1 = -b1;

params(1) = a1;
params(2) = 1-a1;
params(3) = b1;
params(4) = -b1;

%Deterministic steady-state
nxss=[Bh1;Bf1;Sh1;Sf1;zeros(7,1)]; %for states
nyss = zeros(11,1); %for controls
% Make eta
eta=[zeros(4,6);ty 0 0 0 0 0;0 ty 0 0 0 0;0 0 td 0 0 0;0 0 0 td 0 0;0 0 0 0 mtq 0;0 0 0 0 0 mtq;0 0 0 0 0 0];

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices of derivatives H1,..,Hk .
%The external algorithm is Levintal's function solve_dsge.m
algo='gensylv'; %algorithm to compute first-order derivatives
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,2,algo);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0=nxss; % start at the steady state
x0(eps_ind) = 1; % impose the model of interest
derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon follows a unit-root process.

%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
xt = dr_ht(derivs1,nxss,2,(x0-nxss));

%Return h1(0,1,1)-x
eq1 = xt(1)-b1; %first state is home bonds
eq2 = xt(3) - a1; %third state is home equities

y = [eq1;eq2];


