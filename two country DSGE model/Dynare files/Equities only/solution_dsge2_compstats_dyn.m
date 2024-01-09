%--------------------------------------------------------------------------
% Two-country DSGE model with bonds and equities: Comparative statics
% This script replicates (part of) Figure 1 in the paper.
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%--------------------------------------------------------------------------

clear;

addpath('C:\dynare\5.2\matlab');
load('my_dsge2_model.mat');

eps_ind = 3; %index of perturbation variable

%--------------------------------------------------------------------------
%Comparative statics: std. redistributive shocks
%--------------------------------------------------------------------------

mtd = 0.025:0.005:0.1;
l_msd = length(mtd);
sol_msd = zeros(l_msd,1);

tic
for t=1:l_msd
    
    M_.params(10) = mtd(t);
    % this function computes the SSS of equities
    sol_msd(t) = compute_sss_dsge2_dyn(M_,options_,oo_,eps_ind,0.5);
    
    
end
toc

%--------------------------------------------------------------------------
%Comparative statics: std. preference shocks
%--------------------------------------------------------------------------
M_.params(10)=td; % set td back to its calibrated value

mtq = 0:0.0005:0.01;
l_mtq = length(mtq);
sol_mtq = zeros(l_mtq,1);

tic
for t=1:l_mtq
    
    M_.params(14) =mtq(t);
    sol_mtq(t) = compute_sss_dsge2_dyn(M_,options_,oo_,eps_ind,0.5);
    
    
end
toc

% Plot the results
figure;
subplot(1,2,1);
plot(mtd,sol_msd(:,1),'r');
title('(a) As a function of redistributive risk');
xlabel('\eta_{d}');
subplot(1,2,2);
plot(mtq,sol_mtq(:,1),'r');
title('(b) As a function of preference risk');
xlabel('\eta_{q}');
