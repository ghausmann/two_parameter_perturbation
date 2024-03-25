function y = eval_sss_growth(k_bar,model,params,M,eta,eps_ind,approx)
% Neoclassical growth model
% This function implements the SSS algorithm.
%
% The key input is k_bar (DSS target for capital)
% Other inputs are data file 'model', vector of parameters 'params',
% matrices of moments M, the 'eta' matrix, the index of epsilon in vector
% xt, and the approximation order.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given k_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.

betta = params(1);
alpa = params(3);
delta = params(4);

k1 = k_bar(1); % approximation point for capital, in levels
psi = betta*(alpa/(k1^(1-alpa)) + (1-delta)) - 1; %implied value of psi, given k1

% implied steady-state (in logs)
kss = log(k1);
eyss = log(exp(kss)^alpa);
eiss = log(delta*exp(kss));
css = log(exp(eyss) - exp(eiss));
zss = 0;
epsss = 0;

% collect steady-state values
nxss=[kss;zss;epsss];
nyss=[css;eiss;eyss];
params(8) = psi; %change auxiliary parameter value to the implied one.

%STEP 2: Using the new DSS and other function inputs, call an external algorithm
%to compute the matrices H1,..,Hk of derivatives.
%The external algorithm is Levintal's function solve_dsge.m
algo='gensylv'; %algorithm to compute first-order derivatives
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest epsilon=sigma=1.
x0=nxss; % start at the steady state
% impose the model of interest (eps=1), and make sure epsilon is a state
% constant over time.
x0(eps_ind) = 1;
derivs1.hx(eps_ind,eps_ind) = 1;

%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
x1=dr_ht(derivs1,nxss,approx,(x0-nxss));

%Return h1(0,1,1)- log(k_bar)
y = x1(1) - log(k1);