%--------------------------------------------------------------------------
% Multi-asset DSGE model with bonds and equities: impulse responses.
%
% This script replicates figures 8-10 in Appendix E, from the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;

%Set one of these variables to 1:
irf_income_shock = 1; %for income shock
irf_div_shock = 0; %for dividend shock
irf_pref_shock = 0; %for preference shock

%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
%Load the pre-processing data
load('my_dsge4_model.mat');
eps_ind = 5;

a = M_.params(1);

%--------------------------------------------------------------------------
%Model solution
%--------------------------------------------------------------------------

% First, calibrate the std. of preference shocks to match the observed
% Equity home bias, and solve for SSS external bond position.
mysol_calib = compute_calib_dsge4(a,M_,options_,oo_,eps_ind,[0.01 1])

tq = mysol_calib(1)
a1 = a;
b1 = mysol_calib(2)

M_.params(17) = tq;
%SSS check
% mysol = compute_sss_dsge4(M_,options_,oo_,eps_ind,[a b1])

% Recalculate DSS of auxiliary model
Sh0 = a1;
Sf0 = 1-a1;
Bh0 = b1;
Bf0 = -b1;

M_.params(1) = a1;
M_.params(2) = 1-a1;
M_.params(3) = b1;
M_.params(4) = -b1;

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
y0 = yss; % start at the steady state
y0(eps_ind+2)=1; % evaluate at the model of interest (epsilon=1)

T0 = 10;
T = 50;
Tr = 15;

innovations = zeros(6,(T0 + (T-1)));

if irf_income_shock ==1
    innovations(1,T0+1)=1;
elseif irf_div_shock == 1
    innovations(3,T0+1)=1;
elseif irf_pref_shock == 1
    innovations(5,T0+1)=1;
end

myt =simult_(M_,options_,y0,mdr,innovations',3);
xt = myt(3:13,:); %states
yt = [myt(1:2,:);myt(14:22,:)]; %controls

%Compute time-varying moments of returns
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(6,M_.Sigma_e); %monomials to approximate expectations
P = [betta (gama+gap) kappa d0 pf0 zSh0 zSf0 zBh0 zBf0 P0 Ps0]; %other inputs

my_moments =zeros(14,(T0 + (T-1)));
for t=1:(T0 + (T-1))
    my_moments(:,t) = (cross_moments_dsge4(P,yss,mdr,myt(:,t),epsi_nodes,weight_nodes,3));
end

m_Rsh = my_moments(1,2:end-1);
m_Rsf = my_moments(2,2:end-1);
m_Rbh = my_moments(3,2:end-1);
m_Rbf = my_moments(4,2:end-1);
m_Rfull = my_moments(5,2:end-1);

var_Rsh = my_moments(7,2:end-1);
var_Rsf = my_moments(8,2:end-1);
var_Rbh = my_moments(9,2:end-1);
var_Rbf = my_moments(10,2:end-1);
var_Rfull = my_moments(11,2:end-1);

corr_R_full_non_fin = my_moments(14,2:end-1);

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
%preference shock
qt = exp(xt(6,2:end-1));

c_ratio = (ct./c_st);
rer = (pit./pi_st);
ca_y = 100*(zSht.*dSh + zSft.*dSh + zBht.*dBh + zBft.*dBf);
tb_y = 100*(eyt - pit.*ct);

cif = 100*(-zSht.*dSh - zBft.*dBf);
cod = 100*(zSft.*dSf + zBht.*dBh);

Rsh = ( zSht + dt.*eyt ).*((zSht_).^-1);
Rsf = ( zSft + pft.*d_st.*ey_st ).*((zSft_).^-1);
Rs_dif = log(Rsh) -log(Rsf);

Rbh = (piut./zBht_);
Rbf = (piu_st./zBft_);
Rb_dif = log(Rbh) -log(Rbf);

Rsb_dif = log(Rsh) -log(Rbh);

Wealth_h = ( Sht.*zSht + Sft.*zSft + Bht.*zBht + Bft.*zBft); %total
Wealth_h_ = ( Sht_.*zSht_ + Sft_.*zSft_ + Bht_.*zBht_ + Bft_.*zBft_); %total
dWealth_h = Wealth_h - Wealth_h_;
Wealth_f = ( (1-Sht).*zSht + (1-Sft).*zSft + (-Bht).*zBht + (-Bft).*zBft); %total
eq_Wealth = ( Sht.*zSht + Sft.*zSft ); %equity wealth
wealth_ratio = Wealth_h./Wealth_f;
share_foreign_world = zSft./(zSht+zSft);
eq_share_foreign= (Sft.*zSft)./eq_Wealth;
eq_home_bias = (1 - eq_share_foreign./share_foreign_world);
%external gross positions
ext_assets = (Sft.*zSft + Bht.*zBht);
ext_liabilities = ((1-Sht).*zSht + (-Bft).*zBft);
nfa = ext_assets-ext_liabilities;
ext_eq_assets = (Sft.*zSft );
ext_eq_liabilities = ((1-Sht).*zSht );
share_eq_ext_assets = ext_eq_assets./ext_assets;

%Plot the results
time = 0:(Tr);
figure;

s1 = subplot(5,3,1);
if irf_income_shock ==1
    plot(time,eyt(T0:T0+Tr),'LineWidth',2);
    title('Home Output');
elseif irf_div_shock == 1
    plot(time,dt(T0:T0+Tr),'LineWidth',2);
    title('Home dividend share');
elseif irf_pref_shock == 1
    plot(time,qt(T0:T0+Tr),'LineWidth',2);
    title('Home preference');
end

s2 = subplot(5,3,2);
plot(time,rer(T0:T0+Tr),'LineWidth',2);
title('Real exchange rate');

s3 = subplot(5,3,3);
plot(time,ct(T0:T0+Tr),'LineWidth',2);
title('Consumption');

s4 = subplot(5,3,4);
plot(time,c_ratio(T0:T0+Tr),'LineWidth',2);
title('Consumption ratio');

s5 = subplot(5,3,5);
plot(time,tb_y(T0:T0+Tr),'LineWidth',2);
title('Trade balance');

s6 = subplot(5,3,6);
plot(time,ca_y(T0:T0+Tr),'LineWidth',2);
title('Current account');

s7 = subplot(5,3,7);
plot(time,cif(T0:T0+Tr),'LineWidth',2);
title('Capital inflows');

s8 = subplot(5,3,8);
plot(time,cod(T0:T0+Tr),'LineWidth',2);
title('Capital outflows');

s9 = subplot(5,3,9);
plot(time,nfa(T0:T0+Tr),'LineWidth',2);
title('Net foreign assets');

s10 = subplot(5,3,10);
plot(time,wealth_ratio(T0:T0+Tr),'LineWidth',2);
title('Wealth ratio');

s11 = subplot(5,3,11);
plot(time,Rs_dif(T0:T0+Tr),'LineWidth',2);
title('Return stock diff.');

s12 = subplot(5,3,12);
plot(time,Rb_dif(T0:T0+Tr),'LineWidth',2);
title('Return bond diff.');

s13 = subplot(5,3,13);
plot(time,eq_home_bias(T0:T0+Tr),'LineWidth',2);
title('Equity home bias');

s14 = subplot(5,3,14);
plot(time,share_eq_ext_assets(T0:T0+Tr),'LineWidth',2);
title('Equity/(ext. assets) ratio');

s15 = subplot(5,3,15);
plot(time,corr_R_full_non_fin(T0:T0+Tr),'LineWidth',2);
title('Corr (R^w,non fin income)');
