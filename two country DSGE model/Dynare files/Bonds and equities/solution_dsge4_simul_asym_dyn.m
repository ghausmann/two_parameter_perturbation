%--------------------------------------------------------------------------
% Two-country DSGE model with bonds and equities: Stochastic simulations.
% This script replicates the fourth column of Table 2 in the paper.
% Countries are asymmetric in their degree of risk aversion.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;
rng(0); %Fix seed for pseudorandom number generator.
%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.2\matlab');
%Load the pre-processing data
load('my_dsge4_model.mat');

eps_ind = 5; %index of perturbation variable

%Make countries heterogeneous
gap = -0.25;
M_.params(31) = gap;

tq = 0.0067 -3.4135e-05; %calibrated std. of preference shocks
M_.params(17) = tq;

a = M_.params(1);
%--------------------------------------------------------------------------
%Model solution
%--------------------------------------------------------------------------

%Compute SSS of bonds and equities when countries are heterogeneous
mysol = compute_sss_dsge4_asym_dyn(M_,options_,oo_,eps_ind,[a (1-a) b -b])

%Recalculate the new DSS
a1 = mysol(1);
as1 = mysol(2);
b1 = mysol(3);
bs1 = mysol(4);

Sh0 = a1;
Sf0 = as1;
Bh0 = b1;
Bf0 = bs1;

M_.params(1) = a1;
M_.params(2) = as1;
M_.params(3) = b1;
M_.params(4) = bs1;

% Solve for the DSS of the terms of trade with a non-linear solver (see
% Appendix C).
myfun = @(z)dsge_ss_y4(z,[phi alpa md betta],[a1 as1],[b1 bs1]);
pf0 = fzero(myfun,1);

P0 = ( alpa + (1-alpa)*pf0^(1-phi) )^(1/(1-phi));
Ps0 = ( (1-alpa) + alpa*pf0^(1-phi) )^(1/(1-phi));
C0 = (1/P0)*( (1-md) + a1*md + as1*pf0*md + (1-betta)*(b1*P0 + bs1*Ps0) );
Cs0 = (1/Ps0)*( (1-md)*pf0 + (1-a1)*md + (1-as1)*pf0*md  -(1-betta)*(b1*P0 + bs1*Ps0) );

zSf0 = pf0*(betta/(1-betta))*md;
zBh0 = betta*P0;
zBf0 = betta*Ps0;

M_.params(20) = C0;
M_.params(21) = Cs0;
M_.params(22) = P0;
M_.params(23) = Ps0;
M_.params(24) = pf0;
M_.params(26) = zSf0;
M_.params(27) = zBh0;
M_.params(28) = zBf0;

oo_.steady_state(3) = Bh0;
oo_.steady_state(4) = Bf0;
oo_.steady_state(5) = Sh0;
oo_.steady_state(6) = Sf0;

%Compute the perturbation solution:
%The external algorithm is Dynare matlab function resol.m
[mdr, ~, ~, ~] = resol(0, M_, options_, oo_);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
yss = mdr.ys;
x0 = yss; % start at the steady state
x0(eps_ind+2)=1; % evaluate at the model of interest (epsilon=1)

T0 = 10000;
T = 100000;
%draw pseudo-random innovations
innovations = mvnrnd(zeros(6,1),M_.Sigma_e,(T0 + (T-1)))'; % shocks from period 2 to T
%use the Dynare matlab function simult_.m to simulate the economy:
myt =simult_(M_,options_,x0,mdr,innovations',3);
myt = myt(:,T0+1:end);
xt = myt(3:13,T0+1:end); %states
yt = [myt(1:2,T0+1:end);myt(14:22,T0+1:end)]; %controls

%lags, currents, and leads of variables

%Home bond position
Bht_ = xt(1,1:end-2);
Bht = xt(1,2:end-1);
Bhtp = xt(1,3:end);
%Foreign bond position
Bft_ = xt(2,1:end-2);
Bft = xt(2,2:end-1);
Bftp = xt(2,3:end);
%In changes
dBh = Bht-Bht_;
dBf = Bft-Bft_;
%Home equity position
Sht_ = xt(3,1:end-2);
Sht = xt(3,2:end-1);
Shtp = xt(3,3:end);
%Foreign equity position
Sft_ = xt(4,1:end-2);
Sft = xt(4,2:end-1);
Sftp = xt(4,3:end);
%In changes
dSh = Sht-Sht_;
dSf = Sft-Sft_;
%consumption HOME
ct_ = C0*exp(yt(10,1:end-2));
ct = C0*exp(yt(10,2:end-1));
ctp = C0*exp(yt(10,3:end));
%consumption FOREIGN
c_st_ = Cs0*exp(yt(11,1:end-2));
c_st = Cs0*exp(yt(11,2:end-1));
c_stp = Cs0*exp(yt(11,3:end));
%Home bond price
zBht_=zBh0*exp(yt(1,1:end-2));
zBht=zBh0*exp(yt(1,2:end-1));
zBhtp=zBh0*exp(yt(1,3:end));
%Foreign bond price
zBft_ = zBf0*exp(yt(2,1:end-2));
zBft = zBf0*exp(yt(2,2:end-1));
zBftp = zBf0*exp(yt(2,3:end));
%Home stock price
zSht_=zSh0*exp(yt(4,1:end-2));
zSht=zSh0*exp(yt(4,2:end-1));
zShtp=zSh0*exp(yt(4,3:end));
%Foreign stock price
zSft_ = zSf0*exp(yt(5,1:end-2));
zSft = zSf0*exp(yt(5,2:end-1));
zSftp = zSf0*exp(yt(5,3:end));
%price index HOME
pit_ = P0*exp(yt(6,1:end-2));
pit = P0*exp(yt(6,2:end-1));
pitp = P0*exp(yt(6,3:end));
%price index FOREIGN
pi_st_ = Ps0*exp(yt(7,1:end-2));
pi_st = Ps0*exp(yt(7,2:end-1));
pi_stp = Ps0*exp(yt(7,3:end));
%non-adjusted price index HOME
piut_ = P0*exp(yt(8,1:end-2));
piut = P0*exp(yt(8,2:end-1));
piutp = P0*exp(yt(8,3:end));
%non-adjusted price index FOREIGN
piu_st_ = Ps0*exp(yt(9,1:end-2));
piu_st = Ps0*exp(yt(9,2:end-1));
piu_stp = Ps0*exp(yt(9,3:end));
%terms of trade
pft_ = pf0*exp(yt(3,1:end-2));
pft = pf0*exp(yt(3,2:end-1));
pftp = pf0*exp(yt(3,3:end));
%output HOME
eyt_ = exp(xt(10,1:end-2));
eyt = exp(xt(10,2:end-1));
eytp = exp(xt(10,3:end));
%output FOREIGN
ey_st_ = exp(xt(11,1:end-2));
ey_st = exp(xt(11,2:end-1));
ey_stp = exp(xt(11,3:end));
%dividend share HOME
dt_ = (exp(d0+xt(8,1:end-2))./(1+exp(d0+xt(8,1:end-2))));
dt = (exp(d0+xt(8,2:end-1))./(1+exp(d0+xt(8,2:end-1))));
dtp = (exp(d0+xt(8,3:end))./(1+exp(d0+xt(8,3:end))));
%dividend share FOREIGN
d_st_ = (exp(d0+xt(9,1:end-2))./(1+exp(d0+xt(9,1:end-2))));
d_st = (exp(d0+xt(9,2:end-1))./(1+exp(d0+xt(9,2:end-1))));
d_stp = (exp(d0+xt(9,3:end))./(1+exp(d0+xt(9,3:end))));
%relative consumptions
cr = log(ct./c_st);
%real exchange rate, lag, and change

er = log(pit./pi_st);
er_ = log(pit_./pi_st_);

der = er-er_;
%current account and trade balance
ca = 100*(zSht.*dSh + zSft.*dSh + zBht.*dBh + zBft.*dBf);
tb = 100*(eyt - pit.*ct);
%capital inflows and outflows
cif = 100*(-zSht.*dSh - zBft.*dBf);
cod = 100*(zSft.*dSf + zBht.*dBh);
%net and gross flows
net_flows = cif-cod; %by construction equals -ca
gross_flows = cif+cod;
%returns and returns differentials
Rsh = ( zSht + dt.*eyt ).*((zSht_).^-1); %home equity
Rsf = ( zSft + pft.*d_st.*ey_st ).*((zSft_).^-1); %foreign equity
Rbh = (piut./zBht_); %home bond
Rbf = (piu_st./zBft_); %foreign bond
Rs_dif = log(Rsh./Rsf);
Rb_dif = log(Rbh./Rbf);
Rsb_dif = log(Rsh./Rbh);
%wealth and portfolio shares
Wealth = ( Sht.*zSht + Sft.*zSft + Bht.*zBht + Bft.*zBft); %total
eq_Wealth = ( Sht.*zSht + Sft.*zSft ); %equity wealth
wealth_share = Wealth./( zSht + zSft );
% eq_share_home =(Sht.*zSht)./eq_Wealth;
% eq_share_ = eq_share_home(1:end-1);
% eq_share1 = eq_share_home(2:end);
share_foreign_world = zSft./(zSht+zSft);
eq_share_foreign= (Sft.*zSft)./eq_Wealth;
eq_home_bias = (1 - eq_share_foreign./share_foreign_world);
%external gross positions
ext_assets = (Sft.*zSft + Bht.*zBht);
ext_liabilities = ((1-Sht).*zSht + (-Bft).*zBft);
ext_eq_assets = (Sft.*zSft );
ext_eq_liabilities = ((1-Sht).*zSht );
%net foreign assets
nfa = 100*(ext_assets-ext_liabilities);
nfat = nfa(2:end);
nfa_ = nfa(1:end-1);
dnfa = nfat -nfa_;

%stds. of consumption and GDP
std_c = std(log(ct));
std_y = std(log(eyt));

%Compute Euler errors
[n_nodes,epsi_nodes,weight_nodes] = Monomials_1(6,M_.Sigma_e); %monomials to approximate expectations
%Other inputs for the Euler errors function.
%NOTE!
P = [betta (gama+gap) kappa d0 pf0 zSh0 zSf0 zBh0 zBf0 P0 Ps0];

my_errors =zeros(1,T);
for t=1:T
    my_errors(t) = log10(euler_errors_dsge4_dyn(P,yss,mdr,myt(:,t),epsi_nodes,weight_nodes,3));
end
errors_stats = [mean(my_errors) median(my_errors) max(my_errors)]

%results
disp('-----------------------------------------------');
disp('Simulated Moments of the two-country DSGE model')
disp('-----------------------------------------------');
disp('Averages:')
m_eq_home_bias = mean(eq_home_bias)
m_ext_eq_assets = mean(ext_eq_assets)
m_ext_assets = mean(ext_assets)
m_ext_liabilities = mean(ext_liabilities)
m_euler_errors = mean(my_errors)
max_euler_errors = max(my_errors)
disp('Standard deviations:')
std_c_std_y = std_c/std_y
std_er_std_y = std(er)/std_y
std_tb = std(tb)
std_ca = std(ca)
std_dnfa = std(dnfa)
std_grossflows = std(gross_flows)
std_home_bias = std(eq_home_bias)
disp('')
disp('Serial correlations:')
autocorr_tb = corr(tb(2:end)',tb(1:end-1)')
autocorr_ca = corr(ca(2:end)',ca(1:end-1)')
autocorr_dnfa = corr(dnfa(2:end)',dnfa(1:end-1)')
%autocorr_dnfa1 = corr(dnfa1(2:end)',dnfa1(1:end-1)')
autocorr_home_bias = corr(eq_home_bias(2:end)',eq_home_bias(1:end-1)')
disp('')
disp('Cross-correlations:')
corr_er_cr = corr(er',cr')
corr_c_cs = corr(ct',c_st')
corr_cif_cod = corr(cif',cod')
corr_ca_y = corr(ca',eyt')
corr_dnfa_y = corr(dnfa',eyt(2:end)')
disp('')
disp('Conditional hedge ratios:')
my_est = fitlm([Rb_dif' Rs_dif'],der');
my_coefs = table2array(my_est.Coefficients(2:3,1))
ratio_bonds = my_coefs(1)
ratio_equities = my_coefs(2)

