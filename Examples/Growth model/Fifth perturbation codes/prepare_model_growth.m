% Neoclassical growth model.
% Auxiliary model (sigma,eps)=(0,0) uses a discount factor different than 
% the one intended for the model of interest (sigma,eps)=(1,1) to obtain
% approximations around a point different than the DSS.
%
% The model is defined by:
%
% Ef(yp,y,xp,x)=0
% y=g(x);
% xp=h(x)+eta*ep
% h(x)=[tilh(x);Phi(x)];
%
% The functions f and Phi are known.
% The functions g and tilh are unknown.

clear

% Step 1: Define the model by symbolic variables.

% Declare parameters
%(psi is the auxiliary parameter)
syms betta gama alpa delta rho_z sz rho_eps psi real

% Declare state variables: log(capital), log(TFP), epsilon (perturbation
% parameter)
% today
syms k z eps real
% tomorrow
syms kp zp epsp real

% Declare controls: log(consumption), log(investment), log(output)
% today
syms c ei ey real
%tomorrow
syms cp eip eyp real

% Collect state variables:
x=[k, z, eps]; % today
xp=[kp, zp, epsp]; % tomorrow

% Collect control variables: 
y=[c, ei, ey]; % today
yp=[cp, eip, eyp]; % tomorrow

% Collect Parameters
symparams=[betta, gama, alpa, delta, rho_z, sz, rho_eps, psi];

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_z*z;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
%We have 1 zero-mean shock
eta=[0;sz;0];

%Model equations

%budget constraint
f1 = exp(ey) - exp(c) - exp(ei);
%output production function
f2 = exp(z)*(exp(k)^alpa) - exp(ey);
%law of motion capital
f3 = exp(ei) - exp(kp) + (1-delta)*exp(k);
%Euler equation:
%Here ((betta/(1 + psi*(1-eps) ) is the overall discount factor. When
%eps=1, it goes back to just betta.
f4 = (betta/(1 + psi*(1-eps) ))*exp(-gama*cp)*( alpa*exp(eyp-kp) + (1-delta)) ...
    - exp(-gama*c);

f=[f1;f2;f3;f4]; % all conditions

% Choose approximation order
approx=3; % third order solution.

% call differentiate_dsge
model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')