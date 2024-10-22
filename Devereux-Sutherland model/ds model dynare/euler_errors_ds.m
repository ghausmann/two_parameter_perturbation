function y = euler_errors_ds(P,yss,mdr,yt,epsi_nodes,weight_nodes,approx)
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% This function computes Comsuption and Portfolio Euler errors
%
%  Copyright (C) 2024 Guillermo Hausmann Guil

%Paremeters
betta = P(1);
gama = P(2);
kappa = P(3);
d0 = P(4);

zht=yt(1);
zft = yt(2);
ct = yt(13);

%tomorrow's variables (innovations approximated with monomials)
%yt1 is a (nx+ny)-by-N matrix, where N is the number of rows of the matrix
%epsi_nodes. Each column of yt1 are tomorrow's variables under a particular
%realization of future innovations.
yt1 = dr_yt(mdr,yss,approx,yt(5:12,1)-yss(5:12,1),epsi_nodes');

eykt1 = yt1(10,:);
eyk_st1 = yt1(11,:);
dt1 = yt1(12,:);

ct1 = yt1(13,:);
Rzh1 = ( (d0+dt1).*exp(eykt1) ).*(((zht)).^-1);
Rzf1 = ( (d0+dt1).*exp(eyk_st1) ).*(((zft)).^-1);

myexp_int_bh = Rzh1.*(ct1.^-gama);
myexp_int_bf = Rzf1.*(ct1.^-gama);

myexp_bh = weight_nodes'*(myexp_int_bh');
myexp_bf = weight_nodes'*(myexp_int_bf');

errors_ah = abs( ((betta*((ct^-kappa))*myexp_bh )^(-1/gama))/(ct) -1  );
errors_port = abs(myexp_bh/myexp_bf - 1);


y = [errors_ah;errors_port];