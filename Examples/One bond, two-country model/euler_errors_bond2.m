function y = euler_errors_bond2(P,nxss,nyss,x_1,epsi_nodes,weight_nodes,derivs,eta,approx)
% This function computes Euler equation errors of the one bond, two-country model.
% It assumes the model of interest.
% Inputs are P (vector of parameters), nxss and nyss (DSS of states and
% controls in the auxiliary model), x_1 (current states), epsi_nodes (nodes
% of monomials), weight_nodes (probabilities of each node), derivs
% (matrices of derivatives of the policy functions), the eta matrix, and
% approx (order or approximation).

betta = P(1);
gama = P(2);
C0h = P(3);
C0f = P(4);
gap = P(5);

st = x_1-nxss;

%Evaluate decision rules for states and controls
mht = dr_ht(derivs,nxss,approx,st);
mgt = dr_gt(derivs,nyss,approx,st);

Rt = mgt(1); %log interest rate
cht = C0h*mgt(2); %consumption (normalized by DSS)
cft = C0f*mgt(3); %consumption (normalized by DSS)

%Update states with future innovations (approximated with monomials)
my_mht = mht + eta*epsi_nodes';

%Evaluate future controls for each node of the monomials
my_gt1 = dr_gt(derivs,nyss,approx,my_mht-nxss);
cht1 = C0h*my_gt1(2,:);
cft1 = C0f*my_gt1(3,:);

%Compute expectations of the Euler equation
myexp_int_h = ( cht1.^-gama  );
myexp_int_f = ( cft1.^-gama  );

myexp_h = weight_nodes'*(myexp_int_h');
myexp_f = weight_nodes'*(myexp_int_f');

%Return absolute value of errors
errors_h = abs( (((betta+gap)*Rt*myexp_h )^(-1/gama))/(cht) -1  );
errors_f = abs( (((betta-gap)*Rt*myexp_f )^(-1/gama))/(cft) -1  );

y = [errors_h;errors_f];

