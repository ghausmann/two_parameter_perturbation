%--------------------------------------------------------------------------
% Multi-asset DSGE model with bonds and equities: Comparative statics
%
% This script replicates (part of) Figure 6 in Appendix E, from the paper:
% "Solving DSGE models with incomplete markets by perturbation"
% by Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;
%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
%Load the pre-processing data
load('my_dsge4_model.mat');

eps_ind = 5; %index of perturbation variable

%--------------------------------------------------------------------------
%Comparative statics: std. redistributive shocks
%--------------------------------------------------------------------------

mtd = 0.025:0.005:0.1;
l_mtd = length(mtd);
sol_mtd = zeros(l_mtd,2);

tic
for t=1:l_mtd
    
    M_.params(13) = mtd(t);
    % this function computes the SSS of bonds and equities
    sol_mtd(t,:) = compute_sss_dsge4(M_,options_,oo_,eps_ind,[1 1]);
    
end
toc

%--------------------------------------------------------------------------
%Comparative statics: std. preference shocks
%--------------------------------------------------------------------------
M_.params(13)=td; % set td back to its calibrated value

mtq = 0:0.0005:0.01;
l_mtq = length(mtq);
sol_mtq = zeros(l_mtq,2);
guess = [0.5 0];

tic
for t=1:l_mtq
    
    M_.params(17) =mtq(t);
    sol_mtq(t,:) = compute_sss_dsge4(M_,options_,oo_,eps_ind,guess);
    
end
toc

% Plot the results
figure;
subplot(1,2,1);
plot(mtd,sol_mtd(:,1),'b');
hold on;
plot(mtd,sol_mtd(:,2),'g--');
legend('S','B');
title('(a) As a function of redistributive risk');
xlabel('\eta_{d}');
subplot(1,2,2);
plot(mtq,sol_mtq(:,1),'b');
hold on;
plot(mtq,sol_mtq(:,2),'g--');
legend('S','B');
title('(b) As a function of preference risk');
xlabel('\eta_{q}');
