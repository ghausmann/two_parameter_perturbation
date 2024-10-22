%--------------------------------------------------------------------------
% Small open economy: Comparative statics
% Auxiliary model is PAC
% One-parameter perturbation
% Compare the output with what you get in"comp_statics_soe_pac_dynare.m" in the folder
% "Dynare codes\Pac"
%--------------------------------------------------------------------------

clear

% Choose parameter values:
R0 = 1.04; %DSS gross interest rate
betta = (1/R0); %discount factor of the auxiliary model
gama = 4; %risk aversion
rho_y = 0.85; %auto-corr income
rho_r = 0; %auto-corr interest rate
ty = 0.02; %conditional std. of income shocks
tr = 0; %conditional std. of interest rate shocks
uy_ur_corr = 0; %conditional correlation between shocks

psi1 = 1.0000e-05;
psi2 = -1.0000e-04; %controls discount factor of the model of interest 

%Target for sss of NFA
b_bar = 0;
C0 = 1 + (1-(1/R0))*b_bar;

% Make a vector of parameter values:
params = [betta gama rho_y rho_r ty tr R0 b_bar C0 psi1 psi2 uy_ur_corr];


%--------------------------------------------------------------------------
% Comparative statics: income risk
%--------------------------------------------------------------------------
msy = 0.005:0.0025:0.04;
lsy = length(msy);

bsss_sy = zeros(lsy,1);
guess_b = 0;

for t=1:lsy
    
    params_t = params;
    params_t(5) = msy(t);
    my_sss_system = @(x)eval_sss_soe_pac_p1_only_y(x,params_t);
    bsss_sy(t) = fsolve(my_sss_system,guess_b,optimset('TolFun',1e-15));
    guessb = bsss_sy(t);

    
end

%--------------------------------------------------------------------------
% Comparative statics: risk-aversion
%--------------------------------------------------------------------------

mgama = 1:0.1:8;
lg = length(mgama);

bsss_g = zeros(lg,1);
guess_b = 0;

for t=1:lg
    
    params_t = params;
    params_t(2) = mgama(t);
    params(2) = mgama(t);
    my_sss_system = @(x)eval_sss_soe_pac_p1_only_y(x,params_t);
    bsss_g(t) = fsolve(my_sss_system,guess_b,optimset('TolFun',1e-15));
    guessb = bsss_g(t);
    
end

figure;
subplot(1,2,1);
plot(msy,bsss_sy,'g');title('income risk');
subplot(1,2,2);
plot(mgama,bsss_g,'g');title('risk aversion');





