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

%Add Dynare to the search path
addpath('C:\dynare\5.2\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%Use a second-order solution for this exercise
options_.order=2;
options_.k_order_solver=0;

%--------------------------------------------------------------------------
% Comparative statics: PAC location parameter
%--------------------------------------------------------------------------

v_a1 = 0:0.1:1.5;
l_v_a1 = length(v_a1);
pert_sum = zeros(2,l_v_a1);
my_zeros = zeros(1,l_v_a1);

for t=1:l_v_a1
    
    M_.params(1) = v_a1(t); % update PAC location parameter
    %compute the second-order correction factor for Home and
    %Foreign equity
    pert_sum(:,t) = eval2_sss_ds(v_a1(t),M_, options_, oo_,eps_ind);
           
end

%--------------------------------------------------------------------------
% Perturbation solution at SSS
%--------------------------------------------------------------------------
%First, solve for the SSS
a1 = compute_sss_ds(M_,options_,oo_,eps_ind,0)  
%Recalculate DSS of auxiliary model
ah1 = a1;
af1 = -a1;
M_1 = M_;
M_1.params(1) = a1;
%Deterministic steady state
yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;
%Compute derivatives of policy functions
[mdr1, ~, ~, ~] = resol(0, M_1, options_, oo_);

yss = mdr1.ys;
y0 = yss; % start at the steady state
y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)

T0 = 1;
T = 6;
Tr = 5;
%draw pseudo-random innovations
innovations = zeros(4,(T0 + (T-1))); % shocks from period 2 to T
%use the Dynare matlab function simult_.m to simulate the economy:
myt1 =simult_(M_1,options_,y0,mdr1,innovations',2);
aht1 = myt1(5,:); %Home asset
aft1 = myt1(6,:); %Foreign asset

% %--------------------------------------------------------------------------
% % Perturbation solution at a=0
% %--------------------------------------------------------------------------
a2 = 0; 
%Recalculate DSS of auxiliary model
ah2 = a2;
af2 = -a2;
M_2 = M_;
M_2.params(1) = a2;
yss = [zh0;zh0;1;1;ah2;af2;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

[mdr2, ~, ~, ~] = resol(0, M_2, options_, oo_);

yss = mdr2.ys;
y0 = yss; % start at the steady state
y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)

%use the Dynare matlab function simult_.m to simulate the economy:
myt2 =simult_(M_2,options_,y0,mdr2,innovations',2);
aht2 = myt2(5,:); %Home asset
aft2 = myt2(6,:); %Foreign asset
% 
% %--------------------------------------------------------------------------
% % Perturbation solution at a=1.5
% %--------------------------------------------------------------------------
a3 = 1.5; 
%Recalculate DSS of auxiliary model
ah3 = a3;
af3 = -a3;
M_3 = M_;
M_3.params(1) = a3;
yss = [zh0;zh0;1;1;ah3;af3;0;0;0;0;0;0;1;1];
oo_.steady_state = yss;

[mdr3, ~, ~, ~] = resol(0, M_3, options_, oo_);

yss = mdr3.ys;
y0 = yss; % start at the steady state
y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)

%use the Dynare matlab function simult_.m to simulate the economy:
myt3 =simult_(M_3,options_,y0,mdr3,innovations',2);
aht3 = myt3(5,:); %Home asset
aft3 = myt3(6,:); %Foreign asset

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
