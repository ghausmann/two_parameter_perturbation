function y = euler_errors_soe(P1,P2,nxss,nyss,derivs,x_1,epsi_nodes,weight_nodes,approx)
%This function computes Euler equation errors of the SOE 
%It assumes the model of interest. 
%Inputs are P1 (vector of conditional stds.), P2 (vector of parameters),
%nxss (DSS of states in the auxiliary model), %nyss (DSS of controls in the
%auxiliary model), derivs (matrices of derivatives of the policy
%functions), x_1 (current states), epsi_nodes (nodes of monomials),
%weight_nodes (probabilities of each node), and approx (order or
%approximation).
%
% Copyright (C) 2024 Guillermo Hausmann Guil

msigmas = P1;

betta = P2(1);
gama = P2(2);
R0 = P2(3);
C0 = P2(4);

st = x_1-nxss;

%Evaluate decision rules for states and controls
mht = dr_ht(derivs,nxss,approx,st);
mgt = dr_gt(derivs,nyss,approx,st);

zt = st(3); %log interest rate
ct = mgt; %consumption (normalized by DSS)

%Update states with future innovations (approximated with monomials)
my_mht = mht-nxss + [zeros(1,length(weight_nodes));(msigmas'.*epsi_nodes');zeros(1,length(weight_nodes));];

%Evaluate future controls for each node of the monomials
my_gt1 = dr_gt(derivs,nyss,approx,my_mht);
ct1 = my_gt1;

%Compute expectations of the Euler equation
myexp_int_a = ( (C0^-gama)*( ct1.^-gama  )  );
myexpfa = weight_nodes'*(myexp_int_a');
myexp_a = myexpfa*(R0*exp(zt));

%Return absolute value of errors
y = abs( ((betta*myexp_a )^(-1/gama))/(C0*(ct)) -1  );
    
