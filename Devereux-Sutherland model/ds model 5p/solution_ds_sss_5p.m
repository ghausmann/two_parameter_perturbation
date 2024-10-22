%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: comparative statics.
% Standard version with asset holdings
% This script replicates Figure 1 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: SSS condition')
disp('-----------------------------------------------');

clear;

% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
a = 0; %PAC location parameter
betta = 0.96; %discount factor
gama = 2; %risk-aversion
pssi = 5e-6; %PAC control parameter
d0 = 0.167; %Share of capital income
tk = 0.018; %std. capital income
tl = 0.018; %std. labor income
rho_eps = 1; % persistence of epsilon
kappa = 0.007; %Uzawa parameter
rho_y = 0; % persistence income
rho_d = 0; % persistence capital share 
uk_ul_corr = -0.5; %correlation between capital and labor income shocks 
z0 = betta*d0; %DSS asset price

%collect parameters
params=[a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d];
%make eta
eta=[zeros(2,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];

%Moments
%Compute the cross moments:
n_e=4; %number of shocks.
M=gaussian_moments(n_e); %if the shocks are independent standard-normal you can use this function to get the cross moments.
%Update variance-covariance matrix
Sigma = [1 0 uk_ul_corr 0;0 1 0 uk_ul_corr;uk_ul_corr 0 1 0;0 uk_ul_corr 0 1];
M.M2 = Sigma;

eps_ind = 8; %index of epsilon in state vector
%Choose algorithm for first-order solution
%algo='gensylv';
algo='vectorize';
approx = 2; %order of approximation for dynamics

%--------------------------------------------------------------------------
% Comparative statics: PAC location parameter
%--------------------------------------------------------------------------

v_a1 = 0:0.1:1.5;
l_v_a1 = length(v_a1);
pert_sum = zeros(2,l_v_a1);
my_zeros = zeros(1,l_v_a1);

for t=1:l_v_a1
    
    params(1) = v_a1(t); % update PAC location parameter
    %compute the second-order correction factor for Home and
    %Foreign equity
    pert_sum(:,t) = eval2_sss_ds_5p(v_a1(t),model,params,M,eta,eps_ind);
           
end

%--------------------------------------------------------------------------
% Perturbation solution at SSS
%--------------------------------------------------------------------------

%First, solve for the SSS
a_sss = compute_sss_ds_5p(model,params,M,eps_ind,0)
%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
params1 = params;
params1(1) = a_sss;
nxss1=[ah1;af1;zeros(6,1)]; %for states
nyss = [1;1;z0;z0;1;1]; %for controls

%DSS check
%mycheck1=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss1(:);nxss1(:);params1(:)]))
%Compute second-order solution
derivs1=solve_dsge(model,params1,M,eta,nxss1,nyss,approx,algo);

%Simulation
x1=nxss1; %start at the DSS
x1(eps_ind) = 1; %evaluate at the model of interest (epsilon=1)
derivs1.hx(eps_ind,eps_ind) = 1; %make sure epsilon follows a unit-root process.

T0 = 1;
T = 6;
Tr = 5;
%make innovations
innovations = zeros(4,(T0 + (T-1))); % shocks from period 2 to T

%use Levintal's function simul.m to simulate the economy:
[~,xt1]=simul(x1,innovations,nyss,nxss1,eta,derivs1,approx,0,model);
aht1 = xt1(1,:); %Home asset
aft1 = xt1(2,:); %Foreign asset

%--------------------------------------------------------------------------
% Perturbation solution at a=0
%--------------------------------------------------------------------------

%Impose approximation point
a2 = 0;
%Recalculate DSS of auxiliary model
ah2 = a2;
af2 = -a2;
params2 = params;
params2(1) = a2;
nxss2=[ah2;af2;zeros(6,1)]; %for states

%DSS check
%mycheck2=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss2(:);nxss2(:);params2(:)]))
%Compute second-order solution
derivs2=solve_dsge(model,params2,M,eta,nxss2,nyss,approx,algo);

x2=nxss2; %start at the DSS
x2(eps_ind) = 1; %evaluate at the model of interest (epsilon=1)
derivs2.hx(eps_ind,eps_ind) = 1; %make sure epsilon follows a unit-root process.

%use Levintal's function simul.m to simulate the economy:
[~,xt2]=simul(x2,innovations,nyss,nxss2,eta,derivs2,approx,0,model);
aht2 = xt2(1,:);
aft2 = xt2(2,:);

%--------------------------------------------------------------------------
% Perturbation solution at a=1.5
%--------------------------------------------------------------------------

%Impose approximation point
a3 = 1.5;
%Recalculate DSS of auxiliary model
ah3 = a3;
af3 = -a3;
params3 = params;
params3(1) = a3;
nxss3=[ah3;af3;zeros(6,1)]; %for states

%DSS check
mycheck2=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss3(:);nxss3(:);params3(:)]))
%Compute second-order solution
derivs3=solve_dsge(model,params3,M,eta,nxss3,nyss,approx,algo);

x3=nxss3; %start at the steady state
x3(eps_ind) = 1; %evaluate at the model of interest (epsilon=1)
derivs3.hx(eps_ind,eps_ind) = 1; %make sure epsilon follows a unit-root process.

%use Levintal's function simul.m to simulate the economy:
[~,xt3]=simul(x3,innovations,nyss,nxss3,eta,derivs3,approx,0,model);
aht3 = xt3(1,:);
aft3 = xt3(2,:);

%Plot all results
time = 0:(Tr);

figure;
subplot(2,2,1);
plot(v_a1,my_zeros,'k--');
hold on;
plot(v_a1,pert_sum(1,:),'r');
title('(a) Correction factor, home asset');
xlabel('$\bar{a}$','Interpreter','latex');
ylabel('$\frac{1}{2}\left(h_{\varepsilon\varepsilon}+h_{\sigma\sigma}\right)$','Interpreter','latex');

subplot(2,2,2);
plot(v_a1,my_zeros,'k--');
hold on;
plot(v_a1,pert_sum(2,:),'r');
title('(b) Correction factor, foreign asset');
xlabel('$\bar{a}$','Interpreter','latex');
ylabel('$\frac{1}{2}\left(h_{\varepsilon\varepsilon}+h_{\sigma\sigma}\right)$','Interpreter','latex');

subplot(2,2,3);
plot(time,aht1(T0:T0+Tr),'b');
hold on;
plot(time,aht2(T0:T0+Tr),'r--');
plot(time,aht3(T0:T0+Tr),'g-.');
title('(c) Home asset over time');
xlabel('Time');
legend('$\bar{a}=0.747$','$\bar{a}=0$','$\bar{a}=1.5$','Interpreter','latex')

subplot(2,2,4);
plot(time,aft1(T0:T0+Tr),'b');
hold on;
plot(time,aft2(T0:T0+Tr),'r--');
plot(time,aft3(T0:T0+Tr),'g-.');
title('(d) Foreign asset over time');
xlabel('Time');
legend('$\bar{a}=0.747$','$\bar{a}=0$','$\bar{a}=1.5$','Interpreter','latex')



