/*
% The model of this mod file is the main application of the paper 
% "Solving DSGE models with incomplete markets by perturbation",
% by Guillermo Hausmann-Guil.
% Infinite-horizon version of the Coeurdacier-Gourinchas (2016) model.
% Multi-asset DSGE model with bonds and equities.
% Includes income, redistributive, and preference shocks
% Auxiliary model uses PAC

*/

% The code approximates log-deviations of all variables (except bonds and equities, because they can take zero values)
% Hence their deterministic steady-state values in levels enter as parameters in the model equations
% Make sure that the declaration order coincides with the dr order, to avoid a mess!
var zBh, zBf, Bh, Bf, Sh, Sf, eps, q, q_s, d, d_s, ey, ey_s, pf, zSh, zSf, pi, pi_s, piu piu_s, c, c_s;
varexo uy,uy_s,ud,ud_s,uq,uq_s;

parameters a, af, b, bf, betta, gama, pssi, phi, alpa, md, d0, rho_d, td, rho_y, ty, rho_q, tq, rho_eps, kappa, C0, Cs0, P0, Ps0, pf0, zSh0, zSf0, zBh0, zBf0, uy_uys_corr, uy_ud_corr, gap, pssi2;

% Calibration of Table 1 in the paper
a = 0.83; % U.S. equity home bias (this is a target value)
af = 1-a; % by symmetry
b = 0; % to be solved 
bf = -b; % by symmetry
betta = 0.96; % standard value
gama = 2; % standard value 
pssi = 0.0001; % improves accuracy
pssi2 = 0.1; % improves accuracy
phi = 2; % Imbs and Mejean (2015)
alpa = 0.85; % U.S. import share
md = 0.036; % DSS dividend share (FRED data)
d0 = log(md/(1-md));
rho_d = 0.42; % FRED data
td = 0.059; % FRED data
rho_y = 0.51; % FRED data
ty = 0.018; % FRED data
rho_q = 0.46; % within range of rho_d and rho_y
tq = 0.0001; % to be calibrated
rho_eps = 1; % makes perturbation variable constant over time
kappa = 0.007; % match serial correlation of U.S. trade balance
uy_uys_corr = 0.68; % Corsetti, Dedola and Leduc (2008)
uy_ud_corr = 0.12; % FRED data
gap = 0; % gap between risk-aversion parameters across countries (generates asymmetry)

% Deterministic steady-state of the auxiliary model (under perfect symmetry)
% See Appendix E
C0 = 1;
Cs0 = 1;
P0 = 1;
Ps0 = 1;
pf0 = 1;
zSh0 = (betta/(1-betta))*md;
zSf0 = (betta/(1-betta))*md;
zBh0 = betta;
zBf0 = betta;


model;

% exogenous block

% redistributive shocks
d   =  rho_d*(d(-1))   + td*ud;
d_s =  rho_d*(d_s(-1)) + td*ud_s;

% output shocks
ey   =  rho_y*ey(-1)   + ty*uy;
ey_s =  rho_y*ey_s(-1) + ty*uy_s;

% preference shocks
q   =  rho_q*q(-1)   + tq*uq;
q_s =  rho_q*q_s(-1) + tq*uq_s;

% auxiliary perturbation variable
eps = rho_eps*eps(-1);

%endogenous block

%budget constraint HOME
P0*C0*exp(pi+c)     = ( exp(ey)*( 1 - ( exp(d0+d)/(1+exp(d0+d)) )  ) ) 
+ (Sh(-1))*( zSh0*exp(zSh) + exp(ey)*(exp(d0+d)/(1+exp(d0+d)) )) 
+ (Sf(-1))*( zSf0*exp(zSf) + pf0*exp(pf+ey_s)*(exp(d0+d_s)/(1+exp(d0+d_s)) ) ) 
+ P0*(Bh(-1))*exp(piu) + Ps0*(Bf(-1))*exp(piu_s)  
- zSh0*(Sh)*exp(zSh) - zSf0*(Sf)*exp(zSf) 
- zBh0*(Bh)*exp(zBh) - zBf0*(Bf)*exp(zBf);

%budget constraint FOREIGN
%Note that by imposing market-clearing in financial markets, we have:
%Sh_s = 1 - Sh, Sf_s = 1 - Sf, Bh_s = -Bh, Bf_s = -Bf.
Ps0*Cs0*exp(pi_s+c_s) = ( pf0*exp(ey_s+pf)*( 1 - ( exp(d0+d_s)/(1+exp(d0+d_s)) )  ) ) 
+ (1 - Sh(-1))*( zSh0*exp(zSh) + exp(ey)*(exp(d0+d)/(1+exp(d0+d)) )) 
+ (1 - Sf(-1))*( zSf0*exp(zSf) + pf0*exp(pf+ey_s)*(exp(d0+d_s)/(1+exp(d0+d_s)) )) 
+ P0*(-Bh(-1))*exp(piu) + Ps0*(-Bf(-1))*exp(piu_s)
- zSh0*(1-Sh)*exp(zSh) - zSf0*(1-Sf)*exp(zSf)  
- zBh0*(-Bh)*exp(zBh) - zBf0*(-Bf)*exp(zBf);

%Price index HOME
(P0^(1-phi))*exp(pi*(1-phi))   = (    alpa*exp(-(1-phi)*q) + (1-alpa)*(pf0^(1-phi))*exp((pf-q_s)*(1-phi)));
%Price index FOREIGN
(Ps0^(1-phi))*exp(pi_s*(1-phi)) = ((1-alpa)*exp(-(1-phi)*q) +     alpa*(pf0^(1-phi))*exp((pf-q_s)*(1-phi)));

%Non-adjusted price index HOME 
(P0^(1-phi))*exp(piu*(1-phi)) 
= ( alpa + (1-alpa)*(pf0^(1-phi))*exp((pf)*(1-phi)));
%Non-adjusted price index FOREIGN
(Ps0^(1-phi))*exp(piu_s*(1-phi)) 
= ((1-alpa) +  alpa*(pf0^(1-phi))*exp((pf)*(1-phi)));

%Market-clearing HOME good;
alpa*((P0^phi)*C0)*exp(pi*phi + c) + (1-alpa)*((Ps0^phi)*Cs0)*exp(pi_s*phi + c_s) = exp(ey-(phi-1)*q);

%Euler equation Home equity, HOME
betta*(exp(-kappa*c))*( zSh0*exp(zSh(+1)) + exp(ey(+1))*( exp(d0+d(+1))/(1+exp(d0+d(+1))) ) )*exp( -pi(+1)-(gama+gap)*c(+1) ) 
= zSh0*exp(zSh - pi - (gama+gap)*c)*(1 + (1-eps^2)*pssi*(Sh-a));
%Euler equation foreign equity, HOME
betta*(exp(-kappa*c))*( zSf0*exp(zSf(+1)) + pf0*exp(pf(+1)+ey_s(+1))*(exp(d0+d_s(+1))/(1+exp(d0+d_s(+1))) ) )*exp( -pi(+1)-(gama+gap)*c(+1) )
= zSf0*exp(zSf - pi - (gama+gap)*c)*(1 + (1-eps^2)*pssi*(Sf-af));

%Euler equation Home equity, FOREIGN
betta*(exp(-kappa*c_s))*( zSh0*exp(zSh(+1)) + exp(ey(+1))*(exp(d0+d(+1))/(1+exp(d0+d(+1))) ) )*exp( -pi_s(+1)-(gama-gap)*c_s(+1) ) 
= zSh0*exp(zSh - pi_s - (gama-gap)*c_s)*(1 + (1-eps^2)*pssi*((1-Sh)-(1-a)));
%Euler equation foreign equity, FOREIGN
betta*(exp(-kappa*c_s))*( zSf0*exp(zSf(+1)) + pf0*exp(pf(+1)+ey_s(+1))*(exp(d0+d_s(+1))/(1+exp(d0+d_s(+1))) ) )*exp( -pi_s(+1)-(gama-gap)*c_s(+1) )
= zSf0*exp(zSf - pi_s - (gama-gap)*c_s)*(1 + (1-eps^2)*pssi*((1-Sf)-(1-af)));

%Euler equation Home bond, HOME
betta*(exp(-kappa*c))*P0*exp(piu(+1)  -pi(+1) -(gama+gap)*c(+1)) 
= zBh0*exp(zBh-pi-(gama+gap)*c)*(1 + (1-eps^2)*pssi2*(Bh-b));
%Euler equation foreign bond, HOME
betta*(exp(-kappa*c))*Ps0*exp(piu_s(+1)-pi(+1) -(gama+gap)*c(+1)) 
= zBf0*exp(zBf-pi-(gama+gap)*c)*(1 + (1-eps^2)*pssi2*(Bf-bf));

%Euler equation Home bond, FOREIGN
betta*(exp(-kappa*c_s))*P0*exp(piu(+1)-  pi_s(+1)-(gama-gap)*c_s(+1)) 
= zBh0*exp(zBh-pi_s-(gama-gap)*c_s)*(1 + (1-eps^2)*pssi2*(-Bh+b));
%Euler equation foreign bond, FOREIGN
betta*(exp(-kappa*c_s))*Ps0*exp(piu_s(+1)-pi_s(+1)-(gama-gap)*c_s(+1)) 
= zBf0*exp(zBf-pi_s-(gama-gap)*c_s)*(1 + (1-eps^2)*pssi2*(-Bf+bf));

end;

initval;

Sh = a;
Sf = af;
Bh = b;
Bf = bf;
d = 0;
d_s = 0;
ey = 0;
ey_s = 0;
q = 0;
q_s = 0;
eps = 0;
zSh = 0;
zSf = 0;
zBh = 0;
zBf = 0;
pi = 0;
pi_s = 0;
piu = 0;
piu_s = 0;
pf = 0;
c = 0;
c_s = 0;
uy = 0;
uy_s = 0;
ud = 0;
ud_s = 0;
uq = 0;
uq_s = 0;

end;

shocks;
var uy; stderr 1;
var uy_s; stderr 1;
var ud; stderr 1;
var ud_s; stderr 1;
var uq; stderr 1;
var uq_s; stderr 1;
corr uy,uy_s = uy_uys_corr;
corr uy,ud = uy_ud_corr;
corr uy_s,ud_s = uy_ud_corr;
end;

steady; 
check;
stoch_simul(order=3,periods=0,irf=0);
