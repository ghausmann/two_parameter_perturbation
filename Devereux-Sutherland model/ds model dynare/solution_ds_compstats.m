%--------------------------------------------------------------------------
% Devereux-Sutherland (2011) model: Comparative statics
% Standard version with asset holdings
% This script replicates Figure 2 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: comparative statics')
disp('-----------------------------------------------');

clear;

%Add Dynare to the search path
addpath('C:\dynare\5.2\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%--------------------------------------------------------------------------
%Comparative statics: capital income share
%--------------------------------------------------------------------------

v_d = 0.1:0.025:1;
l_v_d = length(v_d);
sol_d = zeros(1,l_v_d);

for t=1:l_v_d
    
    M_d = M_;
    M_d.params(5) = v_d(t); % update income share
    % this function computes the SSS of nominal holdings of the home asset
    sol_d(t) = compute_sss_ds(M_d,options_,oo_,eps_ind,0);  
    
end

%Theoretical solution in DS paper
sol_ds_d = -0.5*(1 + ((1-v_d)./v_d)*uk_ul_corr*(tl/tk) );
%compare results
my_comp_d = [sol_d;sol_ds_d;sol_d-sol_ds_d];

%--------------------------------------------------------------------------
%Comparative statics: correlation between capital and income shocks
%--------------------------------------------------------------------------
v_c = -1:0.1:1;
l_v_c = length(v_c);
sol_c = zeros(1,l_v_c);

for t=1:l_v_c
    
    M_c = M_;
    M_c.params(12) = v_c(t); %update correlation
    uk_ul_corr_c = v_c(t);
    Sigma_c = [1 0 uk_ul_corr_c 0;0 1 0 uk_ul_corr_c;uk_ul_corr_c 0 1 0;0 uk_ul_corr_c 0 1];
    M_c.Sigma_e = Sigma_c;
    M_c.Correlation_matrix = Sigma_c;
    sol_c(t) = compute_sss_ds(M_c,options_,oo_,eps_ind,0);
        
end

%Theoretical solution in DS paper
sol_ds_c = -0.5*(1 + ((1-d0)/d0)*v_c*(tl/tk) );
%compare results
my_comp_c = [sol_c;sol_ds_c;sol_c-sol_ds_c];

%Make figures
figure;
subplot(1,2,1);
plot(v_d,sol_d,'LineWidth',2);
hold on;
plot(v_d,sol_ds_d,'r--','LineWidth',2);
title('(a) As a function of $\overline{\delta}$','Interpreter','latex');
xlabel('$\overline{\delta}$','Interpreter','latex');
legend('SSS portfolio','Zero-order portfolio');

subplot(1,2,2);
plot(v_c,sol_c,'LineWidth',2);
hold on;
plot(v_c,sol_ds_c,'r--','LineWidth',2);
title('(b) As a function of $\mu_{KL}$','Interpreter','latex');
xlabel('$\mu_{KL}$','Interpreter','latex');
legend('SSS portfolio','Zero-order portfolio');


