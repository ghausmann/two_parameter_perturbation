function y = euler_errors_ds_5p(P,nnxss,nnyss,x_1,epsi_nodes,weight_nodes,mdr,eta,approx)
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% This function computes Comsuption and Portfolio Euler errors
% It uses monomials to calculate expectations efficiently.

%Inputs are P (vector of parameters), nnxss (DSS of states),  nnyss (DSS of
%controls), x_1 (current value of states), epsi_nodes (nodes of monomials),
%weight_nodes (probabilities of each node), mdr (derivatives of the policy
%functions), the eta matrix, and approx (order or approximation).

%Paremeters
betta = P(1);
gama = P(2);
kappa = P(3);
d0 = P(4);

st = x_1-nnxss; %states in deviations from DSS

mht = dr_ht(mdr,nnxss,approx,st); %Computes tomorrow's states (before innovations)
mgt = dr_gt(mdr,nnyss,approx,st); %Computes today's controls

ct = mgt(1);
zht=mgt(3);
zft = mgt(4);

%tomorrow's variables (innovations approximated with monomials): my_mht is
%a nx-by-N matrix, where N is the number of rows of the matrix epsi_nodes.
%Each column of my_mht are therefore tomorrow's states under a particular
%realization of future innovations, and each row the value of a given
%tomorrow's state for all possible realizations of innovations.
my_mht = mht + eta*epsi_nodes';
%my_gt1 is a ny-by-N matrix, where N is the number of rows of the matrix
%epsi_nodes. Each column of my_mht are therefore tomorrow's controls under a particular
%realization of future innovations, and each row the value of a given
%tomorrow's control for all possible realizations of innovations. 
my_gt1 = dr_gt(mdr,nnyss,approx,my_mht-nnxss);

eykt1 = my_mht(3,:);
eyk_st1 = my_mht(4,:);
dt1 = my_mht(7,:);
ct1 = my_gt1(1,:);

Rzh1 = ( (d0+dt1).*exp(eykt1) ).*(((zht)).^-1);
Rzf1 = ( (d0+dt1).*exp(eyk_st1) ).*(((zft)).^-1);

myexp_int_bh = Rzh1.*(ct1.^-gama);
myexp_int_bf = Rzf1.*(ct1.^-gama);

%Computing expectations reduces to simple weighted sums:
myexp_bh = weight_nodes'*(myexp_int_bh');
myexp_bf = weight_nodes'*(myexp_int_bf');

errors_ah = abs( ((betta*((ct^-kappa))*myexp_bh )^(-1/gama))/(ct) -1  );
errors_port = abs(myexp_bh/myexp_bf - 1);

y = [errors_ah;errors_port];