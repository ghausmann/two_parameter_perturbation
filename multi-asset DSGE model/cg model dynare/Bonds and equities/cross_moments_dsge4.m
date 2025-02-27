function y = cross_moments_dsge4(P,yss,mdr,yt,epsi_nodes,weight_nodes,approx)
%This function computes conditional moments of asset returns.

%Inputs are P (vector of parameters), yss (DSS), mdr (matrices of
%derivatives of the policy functions), yt (current variables), epsi_nodes
%(nodes of monomials), weight_nodes (probabilities of each node), and
%approx (order or approximation).

%Paremeters
d0 = P(4);
pf0 = P(5);
pSh0 = P(6);
pSf0 = P(7);
pBh0 = P(8);
pBf0 = P(9);
P0 = P(10);
Ps0 = P(11);

%Today's variables
pBht=yt(1);
pBft = yt(2);
pSht=yt(15);
pSft = yt(16);

Bhtp = yt(3);
Bftp = yt(4);
Shtp = yt(5);
Sftp = yt(6);

Wealth = (   (Shtp)*pSh0*exp(pSht) + (Sftp)*pSf0*exp(pSft) + Bhtp*pBh0*exp(pBht) + Bftp*pBf0*exp(pBft) );
prop_sh =((Shtp)*pSh0*exp(pSht))/Wealth;
prop_sf =((Sftp)*pSf0*exp(pSft))/Wealth;
prop_bh =(Bhtp*pBh0*exp(pBht))/Wealth;
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
pSht1 = yt1(15,:);
pSft1 = yt1(16,:);
pit1 = yt1(17,:);
piut1 = yt1(19,:);
piu_st1 = yt1(20,:);
pft1 = yt1(14,:);

Rsh = ( pSh0*exp(pSht1) + (exp(d0+dt1)./(1+exp(d0+dt1))).*exp(eyt1) ).*((pSh0*exp(pSht)).^-1);
Rsf = ( pSf0*exp(pSft1) + pf0*(exp(d0+dst1)./(1+exp(d0+dst1))).*exp(eyst1+pft1) ).*((pSf0*exp(pSft)).^-1); 
Rbh = (P0/pBh0)*exp(piut1-pBht);
Rbf = (Ps0/pBf0)*exp(piu_st1-pBft);
R_full = (prop_sh*Rsh + prop_sf*Rsf + prop_bh*Rbh + prop_bf*Rbf)./(P0*exp(pit1));
non_fin_income = (1- exp(d0+dt1)./(1+exp(d0+dt1))).*(exp(eyt1)./(P0*exp(pit1)));

Rsh2 = Rsh.^2;
Rsf2 = Rsf.^2;
Rbh2 = Rbh.^2;
Rbf2 = Rbf.^2;
Rfull2 = R_full.^2;
R_full_non_fin = R_full.*non_fin_income;
non_fin_income2 = non_fin_income.^2;

%first moments
m_Rsh  = weight_nodes'*(Rsh');
m_Rsf = weight_nodes'*(Rsf');
m_Rbh = weight_nodes'*(Rbh');
m_Rbf = weight_nodes'*(Rbf');
m_Rfull = weight_nodes'*(R_full');
m_non_fin_income = weight_nodes'*(non_fin_income');
m_R_full_non_fin = weight_nodes'*(R_full_non_fin');

%second moments
m2_Rsh  = weight_nodes'*(Rsh2');
m2_Rsf = weight_nodes'*(Rsf2');
m2_Rbh = weight_nodes'*(Rbh2');
m2_Rbf = weight_nodes'*(Rbf2');
m2_Rfull = weight_nodes'*(Rfull2');
m_non_fin_income2 = weight_nodes'*(non_fin_income2');

var_Rsh = m2_Rsh - m_Rsh^2;
var_Rsf = m2_Rsf - m_Rsf^2;
var_Rbh = m2_Rbh - m_Rbh^2;
var_Rbf = m2_Rbf - m_Rbf^2;
var_Rfull = m2_Rfull - m_Rfull^2;
var_non_fin_income = m_non_fin_income2 - m_non_fin_income^2;
cov_R_full_non_fin = m_R_full_non_fin - m_Rfull*m_non_fin_income;
corr_R_full_non_fin = cov_R_full_non_fin/((var_Rfull*var_non_fin_income)^0.5);


y = [m_Rsh m_Rsf m_Rbh m_Rbf m_Rfull m_non_fin_income var_Rsh var_Rsf var_Rbh var_Rbf var_Rfull var_non_fin_income cov_R_full_non_fin corr_R_full_non_fin]';