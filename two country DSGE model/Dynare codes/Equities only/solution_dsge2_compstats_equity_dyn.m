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
l_mtd = length(mtd);
sol_mtd = zeros(l_mtd,1);

for t=1:l_mtd
    
    M_.params(10) = mtd(t); % update std.
    % this function computes the SSS of equities
    sol_mtd(t) = compute_sss_dsge2_dyn(M_,options_,oo_,eps_ind,0.5);
        
end

%--------------------------------------------------------------------------
%Comparative statics: std. preference shocks
%--------------------------------------------------------------------------
M_.params(10)=td; % set td back to its calibrated value

mtq = 0:0.0005:0.01;
l_mtq = length(mtq);
sol_mtq = zeros(l_mtq,1);

for t=1:l_mtq
    
    M_.params(14) =mtq(t); % update std.
    % this function computes the SSS of equities
    sol_mtq(t) = compute_sss_dsge2_dyn(M_,options_,oo_,eps_ind,0.5);
        
end

% Plot the results
figure;
subplot(1,2,1);
plot(mtd,sol_mtd(:,1),'r');
title('(a) As a function of redistributive risk');
xlabel('\eta_{d}');
subplot(1,2,2);
plot(mtq,sol_mtq(:,1),'r');
title('(b) As a function of preference risk');
xlabel('\eta_{q}');
