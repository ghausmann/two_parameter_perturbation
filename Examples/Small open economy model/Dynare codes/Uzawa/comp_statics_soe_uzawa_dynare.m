%---------------------------------------------------------------------------
% Small open economy: Comparative statics
% Auxiliary model is Uzawa
% Remember to run first the script pre_processing_soe_uzawa_basic.m
%--------------------------------------------------------------------------

clear
%Add Dynare's matlab folder to the search path
addpath('C:\dynare\5.5\matlab');
%Load the pre-processing data
load('my_soe_uzawa_basic.mat');

%baseline values of risk-aversion and std. of income shocks
gamma_baseline = 4;
ty_baseline = 0.02;

%index of the epsilon variable in the vector yt (Dynare's endogenous
%variables)
eps_ind = 4;

%--------------------------------------------------------------------------
% Comparative statics: income risk
%--------------------------------------------------------------------------
mty = 0.005:0.0025:0.04;
lty = length(mty);

bsss_ty = zeros(lty,1);
guess_y = 0;

for t=1:lty
    
    M_.params(6) = mty(t);
    % this function computes the SSS of bonds
    bsss_ty(t) = compute_sss_soe_uzawa_dynare(M_,options_,oo_,eps_ind,guess_y);
    guessb = bsss_ty(t);
    
end

%--------------------------------------------------------------------------
% Comparative statics: risk-aversion
%--------------------------------------------------------------------------
M_.params(6) = ty_baseline;

mgama = 1:0.1:8;
lg = length(mgama);

bsss_g = zeros(lg,1);
guess_g = 0;

for t=1:lg
    
    M_.params(2) = mgama(t);
    bsss_g(t) = compute_sss_soe_uzawa_dynare(M_,options_,oo_,eps_ind,guess_g);
    guessb = bsss_g(t);
    
end

M_.params(2) = gamma_baseline;

% Plot the results
figure;
subplot(1,2,1);
plot(mty,bsss_ty,'m');title('income risk');
subplot(1,2,2);
plot(mgama,bsss_g,'m');title('risk aversion');