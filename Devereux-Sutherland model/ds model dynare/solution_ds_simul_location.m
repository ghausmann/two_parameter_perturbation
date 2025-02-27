%--------------------------------------------------------------------------
% Two-country Devereux-Sutherland (2011) model: stochastic simulation.
% Standard version with asset holdings
% This script replicates Figure 5 in the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

disp('-----------------------------------------------');
disp('Devereux-Sutherland model: SSS matters')
disp('-----------------------------------------------');

clear;

rng(0); %Fix seed for pseudorandom number generator.
%Add Dynare to the search path
addpath('C:\dynare\5.5\matlab');
%Load Dynare's model data
load('my_ds_model.mat');

d0 = 1; %Share of capital income (NOTE: only capital income shocks active)
M_.params(5) = d0;
zh0 = betta*d0;
eps_ind = 3; %index of perturbation variable

%--------------------------------------------------------------------------
% Perturbation solution
%--------------------------------------------------------------------------
%First, solve for the SSS

tic;a_sss = compute_sss_ds(M_,options_,oo_,eps_ind,0);toc

%Vector of approximation points for asset holdings
a_v = [-0.8:0.1:-0.6 -0.59:0.01:-0.41 -0.4:0.1:-0.20];
l_a_v = length(a_v);

%Vectors of moments (to be filled in the coming loop)
errors_cee_v = zeros(l_a_v,1);
errors_pee_v = zeros(l_a_v,1);
my_corr = zeros(l_a_v,1);
my_std = zeros(l_a_v,1);

%draw pseudo-random innovations
T0 = 1;
T = 1000;
n_e=4; % number of shocks.
innovations = my_mvnrnd(zeros(n_e,1),M_.Sigma_e,(T0 + (T-1)))'; % shocks from period 2 to T
P = [betta gama kappa d0];
%Use monomials to discretize future innovations
[n_nodes,epsi_nodes,weight_nodes] = Monomials_2(4,M_.Sigma_e); %monomials to approximate expectations


for n = 1:l_a_v
    
    a_sss = a_v(n); %fix the approximation point
    %Recalculate DSS of auxiliary model
    ah1 = a_sss;
    af1 = -a_sss;
    M_.params(1) = a_sss;
    
    %Deterministic steady-state
    yss = [zh0;zh0;1;1;ah1;af1;0;0;0;0;0;0;1;1];
    oo_.steady_state = yss;
    %Compute derivatives
    [mdr, ~, ~, ~] = resol(0, M_, options_, oo_);
            
    %--------------------------------------------------------------------------
    % Simulation of the model of interest
    %--------------------------------------------------------------------------
    yss = mdr.ys;
    y0 = yss; % start at the steady state
    y0(eps_ind+4)=1; % evaluate at the model of interest (epsilon=1)
    
    %use Dynare's function simult_ to simulate the economy
    myt =simult_(M_,options_,y0,mdr,innovations',3);
    myt = myt(:,T0+1:end);
    xt = myt(5:12,:); %states
    yt = [myt(1:4,:);myt(13:end,:)]; %controls
    ct = yt(5,:);
    cft = yt(6,:);
    
    %Compute second moments
    my_smoments = corrcoef(ct,cft);
    my_corr(n) = my_smoments(1,2);
    my_std(n) = std(xt(1,:));
    
    %Loop to compute Euler errors
    my_errors =zeros(2,T);
    for t=1:T
        my_errors(:,t) = (euler_errors_ds(P,yss,mdr,myt(:,t),epsi_nodes,weight_nodes,3));
    end
        
    %Compute mean Euler errors
    errors_cee_v(n) = (mean(my_errors(1,:)));
    errors_pee_v(n) = (mean(my_errors(2,:)));
    
end
%Plot the results
figure;
subplot(1,2,1);
plot(a_v,errors_cee_v);
hold on;
plot(a_v,errors_pee_v','r--');
xlabel('$\bar{a}$','Interpreter','latex');
title('Euler equation errors');
legend('Consumption eq.','Portfolio eq.')
subplot(1,2,2);
plot(a_v,my_corr);
hold on;
plot(a_v,my_std','r--');
xlabel('$\bar{a}$','Interpreter','latex');
title('Second moments');
legend('Corr(ct,cst)','std(aht)')
