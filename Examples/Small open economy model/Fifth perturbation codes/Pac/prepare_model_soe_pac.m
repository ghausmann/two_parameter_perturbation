% Small open economy model with income and interest-rate shocks.
% Auxiliary model is PAC
%The structure of the script follows Levintal (2017).
%You need to run this code only once.
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
syms betta gama rho_y rho_z rho_eps ty tz R0 b_bar C0 psi1 psi2 real
% Declare state variables: bonds, log output, log interest rate, epsilon
% today
syms b logy logz eps real
% tomorrow
syms bp logyp logzp epsp real

% Declare control variables: consumption
%Here c approximates (Ct/C0), where C0 is the DSS of Ct.
% today
syms c real
% tomorrow
syms cp real

% Collect state variables: bonds, (log) output, (log) interest rate
x=[b,logy,logz,eps]; % today 
xp=[bp,logyp,logzp,epsp]; % tomorrow
% Collect control variables: consumption
y=[c]; % today
yp=[cp]; % tomorrow

% Collect parameters
symparams=[betta, gama, rho_y, rho_z, rho_eps, ty, tz, R0, b_bar, C0, psi1, psi2];

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_y*logy;rho_z*logz;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
eta=[0 0;ty 0;0 tz;0 0];

%budget constraint
f1 =   exp(logy) + b ...
    - C0*(c) - ((R0*exp(logz))^(-1))*bp - 0.5*psi1*(1-eps^2)*(bp-b_bar)^2;
%Euler equation
f2 = (betta*(1+psi2*(eps^2)))*(cp^-gama) ...
    - (c^-gama)*( 1/(R0*exp(logz)) + psi1*(1-eps^2)*(bp-b_bar) );

f=[f1;f2]; % all conditions

% Choose approximation order

approx=5; % fifth order solution.

% call differentiate_dsge

model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')