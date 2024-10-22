% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% Auxiliary model uses PAC
%
% Copyright (C) 2024 Guillermo Hausmann Guil
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

%Parameters
syms a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_d;
%States today
syms ah af eyk eyk_s eyl eyl_s d eps real
%States tomorrow
syms ahp afp eykp eyk_sp eylp eyl_sp dp epsp real
%Controls today
syms zh zf c c_s eY eY_s real
%Controls tomorrow
syms zhp zfp cp c_sp eYp eY_sp real

% Collect Parameters
symparams=[a, betta, gama, pssi, d0, tk, tl, rho_eps, kappa, rho_y, rho_d];
% Collect State variables:
x=[ah, af, eyk, eyk_s, eyl, eyl_s, d, eps]; % today
xp=[ahp, afp, eykp, eyk_sp, eylp, eyl_sp, dp, epsp]; % tomorrow
% Collect Control variables: 
y=[c, c_s, zh, zf, eY, eY_s]; % today
yp=[cp, c_sp, zhp, zfp, eYp eY_sp]; % tomorrow

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_y*eyk;rho_y*eyk_s;rho_y*eyl;rho_y*eyl_s;rho_d*d;rho_eps*eps];
% Define eta as in Schmitt-Grohe and Uribe (2004):
%We have 4 zero-mean shocks
eta=[zeros(2,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];

%Model equations

%Endowment HOME
f1 = eY - ( (d0+d)*exp(eyk) + (1-(d0+d))*exp(eyl) );
%Endowment FOREIGN
f2 = eY_s - ( (d0+d)*exp(eyk_s) + (1-(d0+d))*exp(eyl_s) );

%budget constraint HOME
f3 = eY + (ah)*( (d0+d)*exp(eyk) ) + (af)*( (d0+d)*exp(eyk_s) ) ...
- (ahp)*zh - (afp)*zf - (c); 
%budget constraint FOREIGN
f4 = eY_s + (-ah)*( (d0+d)*exp(eyk) ) + (-af)*( (d0+d)*exp(eyk_s) )... 
- (-ahp)*zh - (-afp)*zf - c_s; 

%Euler equation home asset, HOME
f5 = betta*(c^-kappa)*( (d0+dp)*exp(eykp) )*( cp^-gama )... 
- zh*(c^-gama)*(1 + (1-eps^2)*pssi*(ahp-a));
%Euler equation foreign asset, HOME
f6 = betta*(c^-kappa)*( (d0+dp)*exp(eyk_sp) )*( cp^-gama )... 
- zf*(c^-gama)*(1 + (1-eps^2)*pssi*(afp+a));

%Euler equation home asset, FOREIGN
f7 = betta*(c_s^-kappa)*( (d0+dp)*exp(eykp) )*( c_sp^-gama )... 
- zh*(c_s^-gama)*(1 - (1-eps^2)*pssi*(ahp-a));
%Euler equation foreign asset, FOREIGN
f8 = betta*(c_s^-kappa)*( (d0+dp)*exp(eyk_sp) )*( c_sp^-gama )... 
- zf*(c_s^-gama)*(1 - (1-eps^2)*pssi*(afp+a));

f=[f1;f2;f3;f4;f5;f6;f7;f8]; % all conditions

% Choose approximation order
approx=5; % up to a fifth order solution.

% call differentiate_dsge
model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')