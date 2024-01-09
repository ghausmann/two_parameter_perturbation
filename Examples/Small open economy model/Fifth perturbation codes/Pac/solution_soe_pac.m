%---------------------------------------------------------------------------
% Small open economy: Simulation calibrated with mexican data
% Auxiliary model is PAC
%This script computes Euler equation errors and kernel distributions
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%---------------------------------------------------------------------------

clear
rng(0)
% Add folder 'files' to the search path:
addpath('files');
% Load the model:
load('model')

% Choose parameter values:
R0 = 1.0144; %DSS gross interest rate
betta = (1/R0); %discount factor of the auxiliary model
gama = 2; %risk aversion
rho_y = 0.749; %auto-corr income
rho_z = 0.572; %auto-corr interest rate
ty = 0.02723*((1-rho_y^2)^0.5); %conditional std. of income shocks
tz = 0.01958*((1-rho_z^2)^0.5); %conditional std. of interest rate shocks
uy_ur_corr = -0.62; %conditional correlation between shocks

rho_eps = 0.999999; %auto-corr epsilon (make it 1 after computing derivatives)
psi1 = 0.00001; %auxiliary parameter (controls PAC)
psi2 = 0; %controls discount factor of the model of interest (to be calibrated)

%Target for sss of NFA
b_bar = -0.44;
C0 = 1 + (1-(1/R0))*b_bar;

% Make a vector of parameter values:
params = [betta gama rho_y rho_z rho_eps ty tz R0 b_bar C0 psi1 psi2];
% Make eta
eta=[zeros(1,2);ty 0;0 tz;zeros(1,2)];

%DSS values
bss = b_bar;
lyss = 0;
lzss = 0;
epsss = 0;
css = 1;
nxss0=[lyss;lzss;epsss];
nxss=[b_bar;nxss0];
nyss=css;

%Moments
P1 = [ty tz]; %collect stds
% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.
Sigma = [1  uy_ur_corr;uy_ur_corr 1];
M.M2 = Sigma;

algo='gensylv'; % algorithm for the Sylvester equation
eps_ind = 4; %index of epsilon in state vector
approx0 = 2; %order of approximation for SSS
approx1 = 5; %order of approximation for dynamics

%--------------------------------------------------------------------------
% Calibration and perturbation solution
%--------------------------------------------------------------------------
%Calibrate pssi2, and indirectly the discount factor of the model of
%interest:
psi2_calib = compute_betta_soe_pac(P1,b_bar,model,params,M,nxss0,nyss,eps_ind,approx0,0)
params(12) = psi2_calib;
%Just a check
b1_check = compute_sss_soe_pac(P1,model,params,M,nxss0,nyss,eps_ind,approx0,0);

%Compute the perturbation solution:
%First, verify the DSS
dss_check=double(subs(f,[yp(:);y(:);xp(:);x(:);symparams(:)],[nyss(:);nyss(:);nxss(:);nxss(:);params(:)]))
%The external algorithm is Levintal's function solve_dsge.m
derivs=solve_dsge(model,params,M,eta,nxss,nyss,approx1,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------

x0=nxss; % start at the steady state
%Impose the model of interest (epsilon=1)
derivs.hx(eps_ind,eps_ind) = 1;
x0(eps_ind) = 1;

T0 = 1000;
T = 100000;
%draw pseudo-random innovations
innovations = mvnrnd([0 0],Sigma,(T0 + (T-1)))';
%Use Levintal's function simul.m to simulate the economy:
[yt,xt]=simul(x0,innovations,nyss,nxss,eta,derivs,approx1,0,model);
yt = yt(:,T0+1:end);
xt = xt(:,T0+1:end);
ct = C0*yt;

%Compute and plot kernel density functions for NFA and consumption
pd_obj_b = fitdist(xt(1,:)','Kernel','Kernel','epanechnikov');
XXb=-6:0.01:6;
%XX=0:0.01:70;
pd_b_points = pdf(pd_obj_b,XXb);
figure;plot(XXb,pd_b_points,'r');title('wealth');

%--------------------------------------------------------------------------
% Simulation of the pure PAC model (eps=0)
%--------------------------------------------------------------------------

x0_1=nxss; % start at the steady state

[yt_1,xt_1]=simul(x0_1,innovations,nyss,nxss,eta,derivs,approx1,0,model);
yt_1 = yt_1(:,T0+1:end);
xt_1 = xt_1(:,T0+1:end);
ct_1 = C0*yt_1;

%Compute and plot kernel density functions for NFA and consumption
pd_obj_b_1 = fitdist(xt_1(1,:)','Kernel','Kernel','epanechnikov');
%XXb=-6:0.01:6;
XXb_1=0:0.01:70;
pd_b_points_1 = pdf(pd_obj_b_1,XXb_1);

figure;
subplot(1,2,1);
plot(XXb,pd_b_points,'g--');
title('(a) with \epsilon=1')
subplot(1,2,2);
plot(XXb_1,pd_b_points_1,'g--');
title('(b) with \epsilon=0')

%---------------------------------------------------------------------------
%Compute Euler equation errors:
%---------------------------------------------------------------------------
%inputs for the Euler errors function.
P2 = [betta*(1+psi2_calib) gama R0 C0];
[n_nodes1,epsi_nodes,weight_nodes] = Monomials_2(2,Sigma);

lerrors =zeros(1,T);
for t=1:T
    
    lerrors(t) =  log10(euler_errors_soe(P1,P2,nxss,nyss,derivs,xt(:,t),epsi_nodes,weight_nodes,approx1));
    
end

errors_stats = ([mean(lerrors) median(lerrors) max(lerrors)])
figure;histogram(lerrors);title('errors CEE');

