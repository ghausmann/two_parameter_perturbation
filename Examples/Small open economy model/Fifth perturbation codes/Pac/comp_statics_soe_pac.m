%---------------------------------------------------------------------------
% Small open economy: Comparative statics
% Auxiliary model is PAC
%---------------------------------------------------------------------------
clear
% Add folder 'files' to the search path:
addpath('files'); 
% Load the model:
load('model')

% Choose parameter values:
R0 = 1.04; %DSS gross interest rate
betta = (1/R0); %discount factor of the auxiliary model
gama = 4; %risk aversion
rho_y = 0.85; %auto-corr income
rho_z = 0; %auto-corr interest rate
ty = 0.02; %conditional std. of income shocks
tz = 0; %conditional std. of interest rate shocks

rho_eps = 0.999999; %auto-corr epsilon (make it 1 after computing derivatives)
psi1 = 0.00001; %auxiliary parameter (controls PAC)
psi2 = -1.0000e-04; %controls discount factor of the model of interest 

%DSS of NFA
b_bar = 0; %just some initial value
C0 = 1 + (1-(1/R0))*b_bar;

% Make a vector of parameter values:
params = [betta gama rho_y rho_z rho_eps ty tz R0 b_bar C0 psi1 psi2];
% Make eta
eta=[zeros(1,2);ty 0;0 tz;zeros(1,2)];

% DSS values
bss = b_bar;
lyss = 0;
lzss = 0;
epsss = 0;
css = 1;

nxss0=[lyss;lzss;epsss];
nyss=css;

%Moments
% Compute the cross moments:
n_e=size(eta,2); % number of shocks.
M=gaussian_moments(n_e); % if the shocks are independent standard-normal you can use this function to get the cross moments.

% Choose algorithm for the Sylvester equation
algo='gensylv'; % Kamenik algorithm
eps_ind = 4; %index of epsilon in state vector
approx = 2; %order of approximation for SSS

%Approximate SSS of bonds for different perturbation orders
my_sss_2 = compute_sss_soe_pac([ty 0],model,params,M,nxss0,nyss,eps_ind,2,0)
my_sss_3 = compute_sss_soe_pac([ty 0],model,params,M,nxss0,nyss,eps_ind,3,0)
my_sss_4 = compute_sss_soe_pac([ty 0],model,params,M,nxss0,nyss,eps_ind,4,0)
my_sss_5 = compute_sss_soe_pac([ty 0],model,params,M,nxss0,nyss,eps_ind,5,0)

%--------------------------------------------------------------------------
% Comparative statics: income risk
%--------------------------------------------------------------------------
msy = 0.005:0.0025:0.04;
lsy = length(msy);

bsss_sy = zeros(lsy,1);
guess_b = 0;

for t=1:lsy
    
    bsss_sy(t) = compute_sss_soe_pac([msy(t) 0],model,params,M,nxss0,nyss,eps_ind,approx,guess_b);
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
    
    params(2) = mgama(t);
    bsss_g(t) = compute_sss_soe_pac([ty 0],model,params,M,nxss0,nyss,eps_ind,2,guess_b);
    guessb = bsss_g(t);
    
end

figure;
subplot(1,2,1);
plot(msy,bsss_sy,'g');title('income risk');
subplot(1,2,2);
plot(mgama,bsss_g,'g');title('risk aversion');