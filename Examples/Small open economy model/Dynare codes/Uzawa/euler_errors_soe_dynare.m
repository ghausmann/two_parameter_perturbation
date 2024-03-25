function y = euler_errors_soe_dynare(P2,yss,mdr,yt,epsi_nodes,weight_nodes,approx)
% This function computes Euler equation errors of the SOE 
% It assumes the model of interest. 
%
% Inputs are P2 (vector of parameters), yss (DSS), mdr (matrices of
% derivatives of the policy functions), yt (current variables), epsi_nodes
% (nodes of monomials), weight_nodes (probabilities of each node), and
% approx (order or approximation).
%
% Copyright (C) 2024 Guillermo Hausmann Guil

betta = P2(1);
gama = P2(2);
R0 = P2(3);
C0 = P2(4);

%Today's variables
zt = yt(3);
ct = yt(5);

%tomorrow's variables (innovations approximated with monomials)
%yt1 is a (nx+ny)-by-N matrix, where N is the number of rows of the matrix
%epsi_nodes. Each column of yt1 are tomorrow's variables under a particular
%realization of future innovations.
yt1 = dr_yt(mdr,yss,approx,yt(1:4,1)-yss(1:4,1),epsi_nodes');
ct1 = yt1(5,:);

%Compute expectations of the Euler equation
myexp_int_a = ( (C0^-gama)*( ct1.^-gama  )  );
myexpfa = weight_nodes'*(myexp_int_a');
myexp_a = myexpfa*(R0*exp(zt));

%Return absolute value of errors
y = abs( ((betta*myexp_a )^(-1/gama))/(C0*(ct)) -1  );
    
