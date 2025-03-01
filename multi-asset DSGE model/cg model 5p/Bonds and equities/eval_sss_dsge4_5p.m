function y = eval_sss_dsge4_5p(x,model,params,M,eta,eps_ind)
% Multi-asset DSGE model with bonds and equities.
% Baseline calibration with symmetric countries.
% This function implements the SSS algorithm, using a second-order
% approximation (enough for the SSS of this model).
%
% Inputs are x (vector of SSS values), data file 'model', vector of parameters 'params',
% matrices of moments M, the 'eta' matrix, and eps_ind (the index of the
% epsilon variable).

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of bonds and
%equities.
a1 = x(1);
b1 = x(2);

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


