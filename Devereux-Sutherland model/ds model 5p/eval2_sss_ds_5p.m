function y = eval2_sss_ds_5p(a1,model,params,M,eta,eps_ind)
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% This function computes the second-order correction factor for Home and
% Foreign equity, following identical steps as in the SSS algorithm.
%
% The key input is a1 (SSS value of Home equity). 
% Other inputs are data file 'model', vector of parameters 'params',
% matrices of moments M, the 'eta' matrix, and eps_ind (the index of the
% epsilon variable).


%algorithm to compute first-order solution
algo='vectorize';
%algo='gensylv'; 

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of equities.
ah1 = a1;
af1 = -a1;
params(1) = a1;

betta = params(2);
d0 = params(5);
z = betta*d0;

%Deterministic steady state
nxss=[ah1;af1;zeros(6,1)]; %for states
nyss = [1;1;z;z;1;1]; %for controls

%STEP 2: Call an external algorithm to compute the matrices H1,..,Hk .
%The external algorithm is Levintal's function solve_dsge.m
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,2,algo);

%STEP 3: Evaluate the k-order policy rule for states at the DSS, 
%but imposing the model of interest (epsilon=sigma=1).
x0=nxss; %start at the steady-state
x0(eps_ind) = 1; %impose the model of interest
derivs1.hx(eps_ind,eps_ind) = 1; %make sure epsilon follows a unit-root process.

%Evaluate the policy rule
%uses the function dr_ht.m
xt = dr_ht(derivs1,nxss,2,(x0-nxss));

%Return h1(0,1,1)-x (first and second rows only)
y = [xt(1)-a1;xt(2)+a1];