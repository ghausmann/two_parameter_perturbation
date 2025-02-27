%--------------------------------------------------------------------------
% Multi-asset DSGE model with bonds and equities: Stochastic simulations.
%
% This script replicates Table 4 of Appendix E from the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;
rng(0); %Fix seed for pseudorandom number generator.
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

%Approximation choice? (the active choice is set to one)
choice_ba = 1; %baseline (column 1)
choice_fo = 0; %first-order solution (column 2)
choice_so = 0; %second-order solution (column 3)
choice_rl = 0; %risk-linear solution (column 4)
choice_df = 0; %same psi for bonds and equities (column 5)
%NOTE: to replicate column 6, you have to execute first the script
%prepare_model_dsge4_eps1_5p.m, and then set choice_ba = 1;

% Choose parameter values:
a = 0.83; %Just a initial value
af = 1-a;
b = 0; %Just a initial value
bf = -b;
betta = 0.96; % standard value
gama = 2; % standard value
if choice_df == 1
    pssi = 0.001;
    pssi2 = 0.001;
else
    pssi = 0.0001;
    pssi2 = 0.1;
end
phi = 2; % Imbs and Mejean (2015)
alpa = 0.85; % U.S. import share
md = 0.036; %FRED data
d0 = log(md/(1-md));
rho_d = 0.42; %FRED data
td = 0.059; %FRED data
rho_y = 0.51; %FRED data
ty = 0.018; %FRED data
rho_q = 0.46; % within range of rho_d and rho_y
tq = 0.0067; % comparative statics
rho_eps = 1; % makes perturbation variable constant over time
kappa = 0.007; % match serial correlation of U.S. trade balance
%Deterministic steady-state values
C0 = 1;
Cs0 = 1;
P0 = 1;
Ps0 = 1;
pf0 = 1;
zSh0 = (betta/(1-betta))*md;
zSf0 = (betta/(1-betta))*md;
zBh0 = betta;
zBf0 = betta;
gap = 0; %symmetric case (gap>0 introduces assymetric risk aversion)
uy_uys_corr = 0.68; % Corsetti, Dedeloa and Leduc (2008)
uy_ud_corr = 0.12; %FRED data

params=[a af b bf betta gama pssi phi alpa md d0 rho_y ty rho_d td rho_q tq rho_eps kappa C0 Cs0 P0 Ps0 pf0 zSh0 zSf0 zBh0 zBf0 gap pssi2];

%Moments
% Compute the cross moments:
n_e=6; % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 uy_uys_corr uy_ud_corr 0 0 0;uy_uys_corr 1 0 uy_ud_corr 0 0;uy_ud_corr 0 1 0 0 0;0 uy_ud_corr 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
M.M2 = Sigma;

eps_ind = 11; %index of epsilon in state vector
algo='gensylv';
approx = 3; %approximation order to compute derivatives

%approximation order for higher-order solutions
if choice_so == 1
    approx_simul = 2;
else
    approx_simul = 3;
end
 
%Set this variable to 1 for RER not adjusted for preference shocks
not_adjusted = 0;
%--------------------------------------------------------------------------
% Calibration and perturbation solution
%--------------------------------------------------------------------------

% First, calibrate the std. of preference shocks to match the observed
% Equity home bias, and solve for SSS external bond position.
mysol_calib = compute_calib_dsge4_5p(a,model,params,M,eps_ind,[0.01 1]);
             
tq = mysol_calib(1);
a1 = a;
b1 = mysol_calib(2);
params(17) = tq;
eta=[zeros(4,6);ty 0 0 0 0 0;0 ty 0 0 0 0;0 0 td 0 0 0;0 0 0 td 0 0;0 0 0 0 tq 0;0 0 0 0 0 tq;0 0 0 0 0 0];
% % SSS check
% mysol = compute_sss_dsge4_5p(model,params,M,eps_ind,[a b1]);

% Recalculate DSS of auxiliary model
Sh1 = a1;
Sf1 = 1-a1;
Bh1 = b1;
Bf1 = -b1;

params(1) = a1;
params(2) = 1-a1;
params(3) = b1;
params(4) = -b1;

nxss=[Bh1;Bf1;Sh1;Sf1;zeros(7,1)];
nyss = zeros(11,1);

%mycheck=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx_simul,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1;
derivs1.hx(eps_ind,eps_ind) = 1; % evaluate at the model of interest (epsilon=1)

T0 = 10000;
T = 100000;
%draw pseudo-random innovations
innovations = my_mvnrnd(zeros(6,1),Sigma,(T0 + (T-1)))'; % shocks from period 2 to T

%simulate the economy:
if choice_fo == 1
    %First-order case
    [yt,xt]=simul_linear(x0,innovations,nyss,nxss,eta,derivs1);
    %Set proper jacobians and zero-order terms to compute Euler errors:
    %Here these are just plain first-order derivatives
    derivs_hx_c = derivs1.hx;
    derivs_gx_c = derivs1.gx;
    %zero-order terms are the DSS of states and controls
    nnyss = nyss;
    nnxss = nxss;
elseif choice_rl == 1
    %Risk-linear case
    h = 0.0001; %increment to compute numerical derivatives
    [yt,xt]=simul_linear_risk(x0,innovations,nyss,nxss,eta,derivs1,approx,h);
    %Set proper jacobians and zero-order terms to compute Euler errors:
    %Jacobians of a third-order solution to the model of interest
    derivs_hx_c = my_jac_hx(derivs1,nxss,approx,h);
    derivs_gx_c = my_jac_gx(derivs1,nxss,nyss,approx,h);
    %zero-order terms are the SSS of states and controls
    nnyss = dr_gt(derivs1,nyss,approx,(x0-nxss));
    nnxss = dr_ht(derivs1,nxss,approx,(x0-nxss));
else
    [yt,xt]=simul(x0,innovations,nyss,nxss,eta,derivs1,approx_simul,0,model);
end
yt = yt(:,T0+1:end);
xt = xt(:,T0+1:end);

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
dBh = Bhtp-Bht;
dBf = Bftp-Bft;
%Home equity position 
Sht_ = xt(3,1:end-2);
Sht = xt(3,2:end-1);
Shtp = xt(3,3:end);
%Foreign equity position
Sft_ = xt(4,1:end-2);
Sft = xt(4,2:end-1);
Sftp = xt(4,3:end);
%In changes
dSh = Shtp-Sht;
dSf = Sftp-Sft;
%consumption HOME 
ct_ = C0*exp(yt(1,1:end-2));
ct = C0*exp(yt(1,2:end-1));
ctp = C0*exp(yt(1,3:end));
%consumption FOREIGN
c_st_ = Cs0*exp(yt(2,1:end-2));
c_st = Cs0*exp(yt(2,2:end-1));
c_stp = Cs0*exp(yt(2,3:end));
%Home bond price
zBht_=zBh0*exp(yt(3,1:end-2));
zBht=zBh0*exp(yt(3,2:end-1));
zBhtp=zBh0*exp(yt(3,3:end));
%Foreign bond price
zBft_ = zBf0*exp(yt(4,1:end-2));
zBft = zBf0*exp(yt(4,2:end-1));
zBftp = zBf0*exp(yt(4,3:end));
%Home stock price
zSht_=zSh0*exp(yt(5,1:end-2));
zSht=zSh0*exp(yt(5,2:end-1));
zShtp=zSh0*exp(yt(5,3:end));
%Foreign stock price
zSft_ = zSf0*exp(yt(6,1:end-2));
zSft = zSf0*exp(yt(6,2:end-1));
zSftp = zSf0*exp(yt(6,3:end));
%price index HOME
pit_ = P0*exp(yt(7,1:end-2));
pit = P0*exp(yt(7,2:end-1));
pitp = P0*exp(yt(7,3:end));
%price index FOREIGN
pi_st_ = Ps0*exp(yt(8,1:end-2));
pi_st = Ps0*exp(yt(8,2:end-1));
pi_stp = Ps0*exp(yt(8,3:end));
%non-adjusted price index HOME
piut_ = P0*exp(yt(9,1:end-2));
piut = P0*exp(yt(9,2:end-1));
piutp = P0*exp(yt(9,3:end));
%non-adjusted price index FOREIGN
piu_st_ = Ps0*exp(yt(10,1:end-2));
piu_st = Ps0*exp(yt(10,2:end-1));
piu_stp = Ps0*exp(yt(10,3:end));
%terms of trade
pft_ = pf0*exp(yt(11,1:end-2));
pft = pf0*exp(yt(11,2:end-1));
pftp = pf0*exp(yt(11,3:end));
%output HOME
eyt_ = exp(xt(5,1:end-2));
eyt = exp(xt(5,2:end-1));
eytp = exp(xt(5,3:end));
%output FOREIGN
ey_st_ = exp(xt(6,1:end-2));
ey_st = exp(xt(6,2:end-1));
ey_stp = exp(xt(6,3:end));
%dividend share HOME
dt_ = (exp(d0+xt(7,1:end-2))./(1+exp(d0+xt(7,1:end-2))));  
dt = (exp(d0+xt(7,2:end-1))./(1+exp(d0+xt(7,2:end-1))));
dtp = (exp(d0+xt(7,3:end))./(1+exp(d0+xt(7,3:end))));
%dividend share FOREIGN
d_st_ = (exp(d0+xt(8,1:end-2))./(1+exp(d0+xt(8,1:end-2))));
d_st = (exp(d0+xt(8,2:end-1))./(1+exp(d0+xt(8,2:end-1))));
d_stp = (exp(d0+xt(8,3:end))./(1+exp(d0+xt(8,3:end))));
%relative consumptions
cr = log(ct./c_st);
%real exchange rate, lag, and change
if not_adjusted == 0
    er = log(pit./pi_st);
    er_ = log(pit_./pi_st_);
elseif not_adjusted == 1
    er = log(piut./piu_st);
    er_ = log(piut_./piu_st_);
end
der = er-er_;
%current account and trade balance
ca = 100*(zSht.*dSh + zSft.*dSf + zBht.*dBh + zBft.*dBf);
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
Wealth = ( Shtp.*zSht + Sftp.*zSft + Bhtp.*zBht + Bftp.*zBft); %total
eq_Wealth = ( Shtp.*zSht + Sftp.*zSft ); %equity wealth
wealth_share = Wealth./( zSht + zSft );
share_foreign_world = zSft./(zSht+zSft);
eq_share_foreign= (Sftp.*zSft)./eq_Wealth;
eq_home_bias = (1 - eq_share_foreign./share_foreign_world);
%external gross positions
ext_assets = (Sftp.*zSft + Bhtp.*zBht);
ext_liabilities = ((1-Shtp).*zSht + (-Bftp).*zBft);
ext_eq_assets = (Sftp.*zSft );
ext_eq_liabilities = ((1-Shtp).*zSht );
%net foreign assets
nfa = 100*(ext_assets-ext_liabilities);
nfat = nfa(2:end);
nfa_ = nfa(1:end-1);
dnfa = nfat -nfa_;

%stds. of consumption and GDP
std_c = std(log(ct));
std_y = std(log(eyt));

P = [betta (gama+gap) kappa d0 pf0 zSh0 zSf0 zBh0 zBf0 P0 Ps0];
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(6,Sigma);
ln = length(weight_nodes);
my_errors =zeros(1,T);


tic
for t=1:T
    
    if choice_fo ==1 || choice_rl ==1
        my_errors(t) = log10(euler_errors_dsge4_linear_5p(P,nnxss,nnyss,derivs_hx_c,derivs_gx_c,xt(:,t),epsi_nodes,weight_nodes,eta));
    else
        my_errors(t) = log10(euler_errors_dsge4_5p(P,nxss,nyss,xt(:,t),epsi_nodes,weight_nodes,derivs1,eta,approx_simul));
    end
end
toc
errors_stats = ([mean(my_errors) median(my_errors) max(my_errors)])


%results
disp('-----------------------------------------------');
disp('Simulated Moments of the multi-asset DSGE model')
disp('-----------------------------------------------');
disp('Averages:')
m_eq_home_bias = mean(eq_home_bias)
m_ext_eq_assets = mean(ext_eq_assets)
m_ext_assets = mean(ext_assets)
m_ext_liabilities = mean(ext_liabilities)
m_euler_errors = mean(my_errors)
% max_euler_errors = max(my_errors)
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
my_coefs = mvregress([Rb_dif' Rs_dif'],der')
ratio_bonds = my_coefs(1)
ratio_equities = my_coefs(2)


