function y = euler_errors_ds_5p(P,nnxss,nnyss,x_1,epsi_nodes,weight_nodes,mdr,eta,approx)
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


st = x_1-nnxss;

%keyboard;
mht = dr_ht(mdr,nnxss,approx,st);
mgt = dr_gt(mdr,nnyss,approx,st);

ct = mgt(1);
zht=mgt(3);
zft = mgt(4);

%tomorrow's variables (innovations approximated with monomials)
my_mht = mht + eta*epsi_nodes';
my_gt1 = dr_gt(mdr,nnyss,approx,my_mht-nnxss);

eykt1 = my_mht(3,:);
eyk_st1 = my_mht(4,:);
dt1 = my_mht(7,:);
ct1 = my_gt1(1,:);

Rzh1 = ( (d0+dt1).*exp(eykt1) ).*(((zht)).^-1);
Rzf1 = ( (d0+dt1).*exp(eyk_st1) ).*(((zft)).^-1);

myexp_int_bh = Rzh1.*(ct1.^-gama);
myexp_int_bf = Rzf1.*(ct1.^-gama);

myexp_bh = weight_nodes'*(myexp_int_bh');
myexp_bf = weight_nodes'*(myexp_int_bf');

errors_ah = abs( ((betta*((ct^-kappa))*myexp_bh )^(-1/gama))/(ct) -1  );
errors_port = abs(myexp_bh/myexp_bf - 1);


y = [errors_ah;errors_port];