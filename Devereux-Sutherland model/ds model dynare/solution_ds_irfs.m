%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: Impulse responses.
% Standard version with asset holdings
% This script replicates Figure 3 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: impulse responses')
disp('-----------------------------------------------');

clear;

%Add Dynare to the search path
addpath('C:\dynare\5.2\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

rho_d = 0.5; %update persistence of shocks
M_.params(11) = rho_d;
zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%Uncomment for the second-order solution
% options_.order=2;
% options_.k_order_solver=0;

%color of plots?
my_color = 'b';
%my_color = 'r--';
%my_color = 'g:';
%my_color = 'k';

%--------------------------------------------------------------------------
%Perturbation solution
%--------------------------------------------------------------------------
%First, solve for the SSS
a_sss = compute_sss_ds(M_,options_,oo_,eps_ind,0)

% Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
M_1 = M_;
M_1.params(1) = a_sss;
%Deterministic steady-state
yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

%Compute derivatives of policy rules
[mdr, ~, ~, ~] = resol(0, M_1, options_, oo_);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
yss = mdr.ys;
y0 = yss; % start at the steady state
y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)

x0 = y0(5:12); %states
%SSS of all variables
y1 = dr_yt(mdr,yss,2,x0-yss(5:12,1),zeros(4,1));
%Shock to capital income share (note the need to divide by rho_d here,
%since the shock is pre-determined in Dynare's notation)
y1(12) = 0.01/rho_d;

T0 = 1;
T = 50;
Tr = 10;
time_h = 2;

%make innovations
innovations = zeros(4,(T0 + (T-1))); % shocks from period 2 to T

%use Dynare's function simult_ to simulate the economy
myt =simult_(M_,options_,y1,mdr,innovations',3);
%uncomment this one for a second-order solution
%myt =simult_(M_,options_,y1,mdr,innovations',2);
%uncomment this one for a linear-corrected solution
%[myt,~]=simul_linear_risk_dyn(y1(5:12),innovations,mdr,yss,eps_ind,5,12,0.00001);

xt = myt(5:12,:); %states
yt = [myt(1:4,:);myt(13:end,:)]; %controls

%make variables of interest
aht = xt(1,1:end-1);
aft = xt(2,1:end-1);
zht = yt(1,1:end-1);
zft = yt(2,1:end-1);
wealth = zht.*aht + zft.*aft;

%Plot the results
time = 0:(Tr);

figure; 
subplot(2,2,1);
plot(time,zht(T0:T0+Tr),my_color,'LineWidth',2);
title('(a) Home asset price');
subplot(2,2,2);
plot(time,wealth(T0:T0+Tr),my_color,'LineWidth',2);
axis([0 Tr -0.5 0.5]);
title('(b) Wealth');
subplot(2,2,3);
plot(time,aht(T0:T0+Tr),my_color,'LineWidth',2);
title('(c) Home asset holdings');
subplot(2,2,4);
plot(time,aft(T0:T0+Tr),my_color,'LineWidth',2);
title('(c) Foreign asset holdings');

