% Devereux-Sutherland (2011) model
% Version with real wealth, gross rates, and prices as states
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
syms a betta gama pssi d0 tk tl rho_eps kappa rho_y rho_p gap;

%States today
syms wh W zh zf eyk eyk_s eyl eyl_s d eps real
%States tomorrow
syms whp Wp zhp zfp eykp eyk_sp eylp eyl_sp dp epsp real

%Controls today
syms c c_s Rh Rf ah af zh1 zf1 real
%Controls tomorrow
syms cp c_sp Rhp Rfp ahp afp zh1p zf1p real

% Collect State variables:
x=[wh, W, zh, zf, eyk, eyk_s, eyl, eyl_s, d, eps]; % today
xp=[whp, Wp, zhp, zfp, eykp, eyk_sp, eylp, eyl_sp, dp, epsp]; % tomorrow

% Collect Control variables: 
y=[c, c_s, Rh, Rf, ah, af, zh1, zf1]; % today
yp=[cp, c_sp, Rhp, Rfp, ahp, afp, zh1p, zf1p]; % tomorrow

% Collect Parameters
symparams=[a, betta, gama, pssi, d0, tk, tl, rho_eps, kappa, rho_y, rho_p, gap];

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_y*eyk;rho_y*eyk_s;rho_y*eyl;rho_y*eyl_s;rho_p*d;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
%We have 4 zero-mean shocks
eta=[zeros(4,4);tk 0 0 0;0 tk 0 0;0 0 tl 0;0 0 0 tl;0 0 0 0;0 0 0 0];

%Model equations

%budget constraint HOME
f1 = ( (d0+d)*exp(eyk) + (1-(d0+d))*exp(eyl) )... 
+ W*Rf + wh*(Rh-Rf) - Wp - exp(c); 

%budget constraint FOREIGN
f2 = ( (d0+d)*exp(eyk_s) + (1-(d0+d))*exp(eyl_s) )... 
- W*Rf - wh*(Rh-Rf) + Wp - exp(c_s); 

%Equivalent way 1:

%Euler equation home asset, HOME
f3 = (betta+gap*(1-eps^2))*(exp(-kappa*c))*( (d0+dp)*exp(eykp) )*exp( -gama*cp )... 
- zh1*exp( - gama*c)*(1 + (1-eps^2)*pssi*((whp/zh1)-a));
%Euler equation foreign asset, HOME
f4 = (betta+gap*(1-eps^2))*(exp(-kappa*c))*( (d0+dp)*exp(eyk_sp) )*exp( -gama*cp )... 
- zf1*exp( - gama*c)*(1 + (1-eps^2)*pssi*(((Wp-whp)/zf1)+a));

%Euler equation home asset, FOREIGN
f5 = (betta+gap*(1-eps^2))*(exp(-kappa*c_s))*( (d0+dp)*exp(eykp) )*exp( -gama*c_sp )... 
- zh1*exp( - gama*c_s)*(1 + (1-eps^2)*pssi*(-(whp/zh1)+a));
%Euler equation foreign asset, FOREIGN
f6 = (betta+gap*(1-eps^2))*(exp(-kappa*c_s))*( (d0+dp)*exp(eyk_sp) )*exp( -gama*c_sp )... 
- zf1*exp(  - gama*c_s)*(1 + (1-eps^2)*pssi*( -((Wp-whp)/zf1) -a));

%Equivalent way 2:

% %Euler equation home asset, HOME
% f3 = (betta+gap*(1-eps^2))*(exp(-kappa*c))*( d0*exp(eykp) )*exp( -gama*cp )... 
% - zh1*exp( - gama*c)*(1 + (1-eps^2)*pssi*((whp/zh1)-a));
% %Euler equation home asset, FOREIGN
% f4 = (betta+gap*(1-eps^2))*(exp(-kappa*c_s))*( d0*exp(eykp) )*exp( -gama*c_sp )... 
% - zh1*exp( - gama*c_s)*(1 + (1-eps^2)*pssi*(-(whp/zh1)+a));
% %Portfolio Euler equation, HOME
% f5 = (Rhp/(1 + (1-eps^2)*pssi*((whp/zh1)-a)))*exp( -gama*cp ) - (Rfp/(1 + (1-eps^2)*pssi*(((Wp-whp)/zf1)+a)))*exp( -gama*cp );
% %Portfolio Euler equation, FOREIGN
% f6 = (Rhp/(1 + (1-eps^2)*pssi*(-(whp/zh1)+a)))*exp( -gama*c_sp ) - (Rfp/(1 + (1-eps^2)*pssi*( -((Wp-whp)/zf1) -a)))*exp( -gama*c_sp );


%Home and foreign returns
f7 = Rh - (d0+d)*exp(eyk)/((zh));
f8 = Rf - (d0+d)*exp(eyk_s)/((zf));

%Home and foreign asset holdings
f9 = zh1*ah - whp;
f10 = zf1*af - (Wp-whp);

%Asset prices
f11 = zh1 - zhp;
f12 = zf1 - zfp;

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12]; % all conditions

% Choose approximation order
approx=3; % third order solution.

% call differentiate_dsge
model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')