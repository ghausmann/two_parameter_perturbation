% One bond, two-country model.
% The structure of the script follows Levintal (2017).
% You need to run this code only once.
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

clear;

% Step 1: Define the model by symbolic variables.

% Declare parameters
syms betta gama rho_y rho_eps ty b_bar C0h C0f gap psi1 real
% Declare state variables: bonds, log output, log interest rate, epsilon
% today
syms b logyh logyf eps real
% tomorrow
syms bp logyhp logyfp epsp real

% Declare control variables: 
% today
syms R ch cf real
% tomorrow
syms Rp chp cfp real

% Collect state variables: bonds, (log) output, (log) interest rate
x=[b, logyh, logyf, eps]; % today 
xp=[bp, logyhp, logyfp, epsp]; % tomorrow
% Collect control variables: consumption
y=[R, ch, cf]; % today
yp=[Rp, chp, cfp]; % tomorrow

% Collect parameters
symparams=[betta, gama, rho_y, rho_eps, ty, b_bar, C0h, C0f, gap, psi1];

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_y*logyh;rho_y*logyf;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0;ty 0;0 ty;0 0];

%Equilibrium conditions:

%budget constraint HOME
f1 =   exp(logyh) + b - (1/R)*bp - C0h*(ch);
%budget constraint FOREIGN
f2 =   exp(logyf) - b + (1/R)*bp - C0f*(cf);
%Euler equation HOME
f3 = (betta*(1+gap*(eps^2)))*(chp^-gama) ...
    - (ch^-gama)*(1/R)*( 1 + psi1*(1-eps^2)*(bp-b_bar) );
%Euler equation FOREIGN
f4 = (betta*(1-gap*(eps^2)))*(cfp^-gama) ...
    - (cf^-gama)*(1/R)*( 1 - psi1*(1-eps^2)*(bp-b_bar) );


f=[f1;f2;f3;f4]; % all conditions

% Choose approximation order

approx=3; % fifth order solution.

% call differentiate_dsge

model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')