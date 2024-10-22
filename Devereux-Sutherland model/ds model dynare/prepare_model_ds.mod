/*
% Devereux-Sutherland (2011) model
% Standard version with asset holdings
% Auxiliary model uses PAC
%
% Copyright (C) 2024 Guillermo Hausmann Guil
*/

% Make sure that the declaration order coincides with the dr order, to avoid a mess!
var zh, zf, eY, eY_s, ah, af, eps, eyl, eyl_s, eyk, eyk_s, d, c, c_s;
varexo uk,uk_s,ul,ul_s;

parameters a, betta, gama, pssi, d0, tk, tl, rho_eps, kappa, rho_y, rho_d, uk_ul_corr;

% Calibration of Table 1 in the paper
a = 0; %PAC location parameter
betta = 0.96; %discount factor
gama = 2; %risk-aversion
pssi = 5e-6; %PAC control parameter
d0 = 0.167; %Share of capital income
tk = 0.018; %std. capital income
tl = 0.018; %std. labor income
rho_eps = 1; % persistence of epsilon
kappa = 0.007; %Uzawa parameter
rho_y = 0; % persistence income
rho_d = 0; % persistence capital share 
uk_ul_corr = -0.5; %correlation between capital and labor income shocks 

model;

% exogenous block

% capital output shocks
eyk   =  rho_y*eyk(-1) + tk*uk;
eyk_s =  rho_y*eyk_s(-1) + tk*uk_s;
% labor output shocks
eyl   =  rho_y*eyl(-1) + tl*ul;
eyl_s =  rho_y*eyl_s(-1) + tl*ul_s;

% auxiliary perturbation variable
eps = rho_eps*eps(-1);

%weight shock
d = rho_d*d(-1);

%endogenous block

%Endowments
eY = ( (d0+d)*exp(eyk) + (1-(d0+d))*exp(eyl) );
eY_s = ( (d0+d)*exp(eyk_s) + (1-(d0+d))*exp(eyl_s) );

%budget constraint HOME
c  + ah*zh + af*zf   = 
eY + ah(-1)*( (d0+d)*exp(eyk) ) + af(-1)*( (d0+d)*exp(eyk_s) ); 
%budget constraint FOREIGN
c_s  + (-ah)*zh + (-af)*zf   = 
eY_s + (-ah(-1))*( (d0+d)*exp(eyk) ) + (-af(-1))*( (d0+d)*exp(eyk_s) ); 

%Euler equation Home equity, HOME
betta*(c^-kappa)*( (d0+d(+1))*exp(eyk(+1)) )*(c(+1)^-gama) 
= zh*(c^-gama)*(1 + (1-eps^2)*pssi*(ah-a));
%Euler equation foreign equity, HOME
betta*(c^-kappa)*( (d0+d(+1))*exp(eyk_s(+1)) )*(c(+1)^-gama) 
= zf*(c^-gama)*(1 + (1-eps^2)*pssi*(af+a));

%Euler equation Home equity, FOREIGN
betta*(c_s^-kappa)*( (d0+d(+1))*exp(eyk(+1)) )*(c_s(+1)^-gama) 
= zh*(c_s^-gama)*(1 - (1-eps^2)*pssi*(ah-a));
%Euler equation foreign equity, FOREIGN
betta*(c_s^-kappa)*( (d0+d(+1))*exp(eyk_s(+1)) )*(c_s(+1)^-gama) 
= zf*(c_s^-gama)*(1 - (1-eps^2)*pssi*(af+a));

end;

initval;

ah = a;
af = -a;
eY = 1;
eY_s = 1;
zh = betta*d0;
zf = betta*d0;
eyk = 0;
eyk_s = 0;
eyl = 0;
eyl_s = 0;
d = 0;
eps = 0;
c = 1;
c_s = 1;
uk = 0;
uk_s = 0;
ul = 0;
ul_s = 0;

end;

shocks;
var uk; stderr 1;
var uk_s; stderr 1;
var ul; stderr 1;
var ul_s; stderr 1;
corr uk,ul = uk_ul_corr;
corr uk_s,ul_s = uk_ul_corr;
end;

steady; 
check;

stoch_simul(order=3,periods=0,irf=0);
