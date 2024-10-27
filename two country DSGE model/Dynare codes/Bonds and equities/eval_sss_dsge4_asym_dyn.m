function y = eval_sss_dsge4_asym_dyn(x, M_, options_, oo_, eps_ind)
% Two-country DSGE model with bonds and equities
% General case that allows for asymmetry in risk aversion.
% This function implements the algorithm described in Section 2.4 of the
% paper.
% Inputs are x (vector of SSS values), Dynare's structs, and eps_ind (the
% index of the epsilon variable)
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%STEP 1: Given x, recalculate the new DSS and implied auxiliary parameter
%values, modifying other function inputs when necessary.
a1 = x(1);
as1 = x(2);
b1 = x(3);
bs1 = x(4);

Sh1 = a1;
Sf1 = as1;
Bh1 = b1;
Bf1 = bs1;

betta = M_.params(5);
phi = M_.params(8);
alpa = M_.params(9);
md = M_.params(10);

% Solve for the DSS of the terms of trade with a non-linear solver (see
% Appendix C).
myfun = @(z)dsge_ss_y4(z,[phi alpa md betta],[a1 as1],[b1 bs1]);
pf1 = fzero(myfun,1);

P1 = ( alpa + (1-alpa)*pf1^(1-phi) )^(1/(1-phi));
Ps1 = ( (1-alpa) + alpa*pf1^(1-phi) )^(1/(1-phi));
C1 = (1/P1)*( (1-md) + a1*md + as1*pf1*md + (1-betta)*(b1*P1 + bs1*Ps1) );
Cs1 = (1/Ps1)*( (1-md)*pf1 + (1-a1)*md + (1-as1)*pf1*md  -(1-betta)*(b1*P1 + bs1*Ps1) );

pSf1 = pf1*(betta/(1-betta))*md;
pBh1 = betta*P1;
pBf1 = betta*Ps1;

M_.params(1) = a1;
M_.params(2) = as1;
M_.params(3) = b1;
M_.params(4) = bs1;

M_.params(20) = C1;
M_.params(21) = Cs1;
M_.params(22) = P1;
M_.params(23) = Ps1;
M_.params(24) = pf1;
M_.params(26) = pSf1;
M_.params(27) = pBh1;
M_.params(28) = pBf1;

oo_.steady_state(3) = Bh1;
oo_.steady_state(4) = Bf1;
oo_.steady_state(5) = Sh1;
oo_.steady_state(6) = Sf1;

yss = oo_.steady_state;

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices of derivatives H1,..,Hk .
%The external algorithm is Dynare's MATLAB function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest (epsilon=sigma=1).
x0 = yss(3:13,1);
x0(eps_ind) = 1; %This imposes epsilon=1

%This function evaluates the k-order decision rule h_k.
x1 = dr_yt(mdr,yss,2,x0-yss(3:13,1),zeros(6,1));

eq1 =  x1(3)-b1;
eq2 = x1(4)-bs1;
eq3 = x1(5)-a1;
eq4 = x1(6)-as1;

%keyboard;
y = [eq1;eq2;eq3;eq4];