function y = euler_errors_dsge4_linear_5p(P,hss,gss,derivs_hx_c,derivs_gx_c,x_1,epsi_nodes,weight_nodes,eta)
%Paremeters
betta = P(1);
gama = P(2);
kappa = P(3);
d0 = P(4);
pf0 = P(5);
pSh0 = P(6);
pSf0 = P(7);
pBh0 = P(8);
pBf0 = P(9);
P0 = P(10);
Ps0 = P(11);

st = x_1-hss;

mht = dr_ht_1_risk_linear(derivs_hx_c,hss,st);
mgt = dr_gt_1_risk_linear(derivs_gx_c,gss,st);

ct = mgt(1);
zBht=mgt(3);
zBft = mgt(4);
zSht=mgt(5);
zSft = mgt(6);
pit = mgt(7);


Bhtp = mht(1);
Bftp = mht(2);
Shtp = mht(3);
Sftp = mht(4);

my_mxtp = mht + eta*epsi_nodes';
my_gt1 = dr_gt_1_risk_linear(derivs_gx_c,gss,(my_mxtp-hss));

eyt1=my_mxtp(5,:);
eyst1=my_mxtp(6,:);
dt1=my_mxtp(7,:);
dst1=my_mxtp(8,:);

%keyboard;
ct1 = my_gt1(1,:);
zSht1 = my_gt1(5,:);
zSft1 = my_gt1(6,:);
pit1 = my_gt1(7,:);
piut1 = my_gt1(9,:);
piu_st1 = my_gt1(10,:);
pft1 = my_gt1(11,:);

%keyboard;

Wealth = ( (Shtp)*pSh0*exp(zSht) + (Sftp)*pSf0*exp(zSft) + Bhtp*pBh0*exp(zBht) + Bftp*pBf0*exp(zBft) );
prop_sh =((Shtp)*pSh0*exp(zSht))/Wealth;
prop_sf =((Sftp)*pSf0*exp(zSft))/Wealth;
prop_bh =(Bhtp*pBh0*exp(zBht))/Wealth;
prop_bf = 1-(prop_sh+prop_sf+prop_bh);


Rsh1 = ( pSh0*exp(zSht1) + (exp(d0+dt1)./(1+exp(d0+dt1))).*exp(eyt1) ).*((pSh0*exp(zSht)).^-1);
Rsf1 = ( pSf0*exp(zSft1) + pf0*(exp(d0+dst1)./(1+exp(d0+dst1))).*exp(eyst1+pft1) ).*((pSf0*exp(zSft)).^-1); 
Rbh1 = (P0/pBh0)*exp(piut1-zBht);
Rbf1 = (Ps0/pBf0)*exp(piu_st1-zBft);
R_full = prop_sh*Rsh1 + prop_sf*Rsf1 + prop_bh*Rbh1 + prop_bf*Rbf1;

myexp_int_full = R_full.*exp(pit-pit1).*exp(-gama*ct1);
myexp_full = weight_nodes'*(myexp_int_full');

y = abs( ((betta*(exp(-kappa*ct))*myexp_full )^(-1/gama))/exp(ct) -1  );

