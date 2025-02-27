function y = eval_sss_ds_wealth_5p(x,model,params,M,eta,eps_ind)
% Devereux-Sutherland (2011) model
% Version with real wealth, gross rates, and prices as states
% This function implements the SSS algorithm, using a second-order
% approximation (enough for the SSS of this model).
%
% The key input is the 1-by-2 vector (SSS value of real holdings of home asset, and home asset price). 
% Other inputs are data file 'model', vector of parameters 'params',
% matrices of moments M, the 'eta' matrix, and eps_ind (the index of the
% epsilon variable).


algo='gensylv'; %algorithm to compute first-order derivatives

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of equities.
wh1 = x(1);
z = x(2);

betta = params(2);
d0 = params(5);
gap = z/d0 - betta;
%z = (betta+gap)*d0;
R = 1/(betta+gap);

%params(1) = wh1;
params(1) = wh1/z;
params(12) = gap;

%Deterministic steady-state
nxss=[wh1;0;z;z;zeros(6,1)]; %for states
nyss = [zeros(2,1);R;R;wh1/z;-wh1/z;z;z]; %for controls


%keyboard;
%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices of derivatives H1,..,Hk .
%The external algorithm is Levintal's function solve_dsge.m
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,2,algo);



%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0=nxss; % start at the steady state
x0(eps_ind) = 1; % impose the model of interest
derivs1.hx(eps_ind,eps_ind) = 1; % make sure epsilon follows a unit-root process.

%Compute the Taylor-series policy rule for bonds
%uses the function dr_ht.m
xt = dr_ht(derivs1,nxss,2,(x0-nxss));

%Return h1(0,1,1)-a1
eq1 = xt(1)-wh1;
eq2 = xt(3) - z;
%keyboard;
y = [eq1;eq2];
