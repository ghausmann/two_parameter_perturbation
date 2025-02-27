function y = eval2_sss_ds(a1,M_, options_, oo_,eps_ind)
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% This function computes the second-order correction factor for Home and
% Foreign equity, following identical steps as in the SSS algorithm.
%
% The key input is a1 (SSS value of Home equity). 
% Other inputs are Dynare's structs, and eps_ind (the
% index of the epsilon variable)


%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
%In this application, we only need to recalculate the DSS of equities.
ah1 = a1;
af1 = -a1;
M_.params(1) = a1;

betta = M_.params(2);
d0 = M_.params(5);
zh0 = betta*d0;

%Deterministic steady state
yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

%STEP 2: Call an external algorithm to compute the matrices H1,..,Hk .
%The external algorithm is Dynare's MATLAB function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Evaluate the k-order policy rule for states at the DSS, 
%but imposing the model of interest (epsilon=sigma=1).
x0 = yss(5:12,1); %start at the steady-state
x0(eps_ind) = 1; %impose the model of interest

%Evaluate the policy rule
%uses the function dr_yt.m
y1 = dr_yt(mdr,yss,2,x0-yss(5:12,1),zeros(4,1));

%Return h1(0,1,1)-x (rows for home and foreign assets only)
y = [y1(5) - a1;y1(6) + a1];

