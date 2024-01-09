/*
 * Small open-economy model: calibration for simulations with mexican data
 * Auxiliary model is Uzawa
%
% Copyright (C) 2024 Guillermo Hausmann Guil
*/

%Here c approximates (Ct/C0), where C0 is the DSS of Ct.
%Make sure that the declaration order coincides with the dr order, to avoid a mess!
var b,logy,logz,eps,c;
varexo u_y,u_r;

parameters betta, gama, rho_y, rho_z, rho_eps, ty, tz, R0, b_bar, C0, psi1, psi2, uy_ur_corr;

R0 = 1.0144; %DSS gross interest rate
betta = (1/R0); %discount factor of the auxiliary model
gama = 2; %risk aversion
rho_y = 0.749; %auto-corr income
rho_z = 0.572; %auto-corr interest rate
ty = 0.02723*((1-rho_y^2)^0.5); %conditional std. of income shocks
tz = 0.01958*((1-rho_z^2)^0.5); %conditional std. of interest rate shocks
uy_ur_corr = -0.62; %conditional correlation of shocks

rho_eps = 1;
psi1 = 0.00001;
psi2 = 0; %just a initial value (to be calibrated)

%Target for sss of NFA
b_bar = -0.44;
C0 = 1 + (1-(1/R0))*b_bar;

model;

%Laws of motion for the exogenous states
%Output
logy = rho_y*logy(-1) + ty*u_y;
%Interest rate
logz = rho_z*logz(-1) + tz*u_r;
%epsilon variable
eps = rho_eps*eps(-1);


%budget constraint
0 =  exp(logy) + b(-1) - C0*c - ((R0*exp(logz))^(-1))*b;
%Euler equation
0 = (betta*(1+psi2*(eps^2)))*(c(+1)^-gama)*((c)^-(psi1*(1-eps^2))) 
    - (1/R0)*exp(-logz)*(c^-gama);

end;

initval;
b = b_bar;
logy = 0;
logz = 0;
eps = 0;
c = 1;
u_y = 0;
u_r = 0;
end;

shocks;

var u_y; stderr 1;
var u_r; stderr 1;
corr u_y,u_r = uy_ur_corr;

end;

steady;
check;

stoch_simul(order=3,periods=0,irf=0);
