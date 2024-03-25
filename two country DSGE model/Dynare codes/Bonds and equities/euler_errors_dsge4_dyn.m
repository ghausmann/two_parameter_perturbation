function y = euler_errors_dsge4_dyn(P,yss,mdr,yt,epsi_nodes,weight_nodes,approx)
%This function computes Euler equation errors of the Two-country DSGE
%model, using the wealth Euler equation (4.24) in the paper.
%
%Inputs are P (vector of parameters), yss (DSS), mdr (matrices of
%derivatives of the policy functions), yt (current variables), epsi_nodes
%(nodes of monomials), weight_nodes (probabilities of each node), and
%approx (order or approximation).
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%Paremeters
betta = P(1);
gama = P(2);
theta = P(3);
d0 = P(4);
pf0 = P(5);
pSh0 = P(6);
pSf0 = P(7);
pBh0 = P(8);
pBf0 = P(9);
P0 = P(10);
Ps0 = P(11);

%Today's variables
ct = yt(21);
zBht=yt(1);
zBft = yt(2);
zSht=yt(15);
zSft = yt(16);
pit = yt(17);

Bhtp = yt(3);
Bftp = yt(4);
Shtp = yt(5);
Sftp = yt(6);

Wealth = (   (Shtp)*pSh0*exp(zSht) + (Sftp)*pSf0*exp(zSft) + Bhtp*pBh0*exp(zBht) + Bftp*pBf0*exp(zBft) );
prop_sh =((Shtp)*pSh0*exp(zSht))/Wealth;
prop_sf =((Sftp)*pSf0*exp(zSft))/Wealth;
prop_bh =(Bhtp*pBh0*exp(zBht))/Wealth;
prop_bf = 1-(prop_sh+prop_sf+prop_bh);

%tomorrow's variables (innovations approximated with monomials)
%yt1 is a (nx+ny)-by-N matrix, where N is the number of rows of the matrix
%epsi_nodes. Each column of yt1 are tomorrow's variables under a particular
%realization of future innovations.
yt1 = dr_yt(mdr,yss,approx,yt(3:13,1)-yss(3:13,1),epsi_nodes');

dt1=yt1(10,:);
dst1=yt1(11,:);
eyt1=yt1(12,:);
eyst1=yt1(13,:);
ct1 = yt1(21,:);
pSht1 = yt1(15,:);
pSft1 = yt1(16,:);
pit1 = yt1(17,:);
piut1 = yt1(19,:);
piu_st1 = yt1(20,:);
pft1 = yt1(14,:);

Rsh1 = ( pSh0*exp(pSht1) + (exp(d0+dt1)./(1+exp(d0+dt1))).*exp(eyt1) ).*((pSh0*exp(zSht)).^-1);
Rsf1 = ( pSf0*exp(pSft1) + pf0*(exp(d0+dst1)./(1+exp(d0+dst1))).*exp(eyst1+pft1) ).*((pSf0*exp(zSft)).^-1); 
Rbh1 = (P0/pBh0)*exp(piut1-zBht);
Rbf1 = (Ps0/pBf0)*exp(piu_st1-zBft);
R_full = prop_sh*Rsh1 + prop_sf*Rsf1 + prop_bh*Rbh1 + prop_bf*Rbf1;

%Compute expectations of the Euler equation
myexp_int_a = R_full.*exp(pit-pit1).*exp(-gama*ct1);
myexp_a = weight_nodes'*(myexp_int_a');


y = abs( ((betta*(exp(-theta*ct))*myexp_a )^(-1/gama))/exp(ct) -1  );
