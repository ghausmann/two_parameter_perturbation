%--------------------------------------------------------------------------
% Small open economy: Simulation calibrated with mexican data
% Auxiliary model is PAC
% This script performs a stochastic simulation of the model.
% One-parameter perturbation approach
% Compare the output with what you get in "solution_soe_pac_dynare_testing.m" in the folder
% "Dynare codes\Pac"
%--------------------------------------------------------------------------

clear
rng(0); %fix seed for full replication

% Choose parameter values:
R0 = 1.0144; %DSS gross interest rate
betta = (1/R0); %discount factor of the auxiliary model
gama = 2; %risk aversion
rho_y = 0.749; %auto-corr income
rho_r = 0.572; %auto-corr interest rate
ty = 0.02723*((1-rho_y^2)^0.5); %conditional std. of income shocks
tr = 0.01958*((1-rho_r^2)^0.5); %conditional std. of interest rate shocks
uy_ur_corr = -0.62; %conditional correlation between shocks

psi1 = 0.00001; %auxiliary parameter (controls PAC)
psi2 = 0; %controls discount factor of the model of interest (to be calibrated)

%Target for sss of NFA
b_bar = -0.44;
C0 = 1 + (1-(1/R0))*b_bar;

% Make a vector of parameter values:
params = [betta gama rho_y rho_r ty tr R0 b_bar C0 psi1 psi2 uy_ur_corr];

%Solve for the value of psi2 that lead to the target SSS for bonds
my_calib_system = @(x)calib_sss_soe_pac_p1(x,b_bar,params);
my_psi2_sol = fsolve(my_calib_system,0)
params(11) = my_psi2_sol;

%Check the solution
my_sss_system = @(x)eval_sss_soe_pac_p1(x,params);
my_sss_sol = fsolve(my_sss_system,b_bar)

%Compute the third-order solution
sol_o1 = sol_soe_o1(params)
my_o2_system = @(x)system_soe_o2(x,sol_o1,params);
sol_o2 = fsolve(my_o2_system,zeros(1,7))
my_o3_system = @(x)system_soe_o3(x,sol_o1,sol_o2,params);
sol_o3 = fsolve(my_o3_system,zeros(1,13))

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------

T0 = 100;
T = 1000;
Sigma = [1 uy_ur_corr;uy_ur_corr 1];
%draw pseudo-random innovations
innovations = mvnrnd([0 0],Sigma,(T0 + (T-1)))';

%states
btp = zeros(1,T0+T);
yt = zeros(1,T0+T);
rt = zeros(1,T0+T);

btp(1) = b_bar;

for t=2:T0+T

    rt(t) = rho_r*rt(t-1) + tr*innovations(2,t-1);
    yt(t) = rho_y*yt(t-1) + ty*innovations(1,t-1);
    St = [btp(t-1) rt(t) yt(t)];
    btp(t) = btp_eval(b_bar,St,sol_o1,sol_o2,sol_o3);

end

btp_comp = btp(T0+1:end);

figure;plot(btp_comp, 'r--');


