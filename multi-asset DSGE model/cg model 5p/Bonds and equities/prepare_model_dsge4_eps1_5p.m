% The model generated by this script is the main application of the paper 
% "Solving DSGE Models with Incomplete Markets by Perturbation",
% by Guillermo Hausmann-Guil.
%
% Infinite-horizon version of the Coeurdacier-Gourinchas (2016) model.
% Multi-asset DSGE model with bonds and equities.
% Includes income, redistributive, and preference shocks
% Auxiliary model uses PAC
% NOTE: in this model, the perturbation parameter epsilon is NOT squared,
% but everything else is identical to the script prepare_model_dsge4_5p.m

%
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

%Parameters.
% The code approximates log-deviations of all variables (except for bonds and equities, because they can take zero values)
% Hence their deterministic steady-state values in levels enter as parameters in the model equations
syms a af b bf betta gama pssi phi alpa md d0 rho_d td rho_y ty rho_q tq rho_eps kappa C0 Cs0 P0 Ps0 pf0 zSh0 zSf0 zBh0 zBf0 gap pssi2;

%States today
syms Bh Bf Sh Sf d d_s ey ey_s q q_s eps real
%States tomorrow
syms Bhp Bfp Shp Sfp dp d_sp eyp ey_sp qp q_sp epsp real

%Controls today
syms c c_s zBh zBf zSh zSf pi pi_s piu piu_s pf real
%Controls tomorrow
syms cp c_sp zBhp zBfp zShp zSfp pip pi_sp piup piu_sp pfp real

% Collect State variables:
x=[Bh, Bf, Sh, Sf, ey, ey_s, d, d_s, q, q_s, eps]; % today
xp=[Bhp, Bfp, Shp, Sfp, eyp, ey_sp, dp, d_sp, qp, q_sp, epsp]; % tomorrow

% Collect Control variables: 
y=[c, c_s, zBh, zBf, zSh, zSf, pi, pi_s, piu, piu_s, pf]; % today
yp=[cp, c_sp, zBhp, zBfp, zShp, zSfp, pip, pi_sp, piup, piu_sp, pfp]; % tomorrow

% Collect Parameters
symparams=[a, af, b, bf, betta, gama, pssi, phi, alpa, md, d0, rho_y, ty, rho_d, td, rho_q, tq, rho_eps, kappa, C0, Cs0, P0, Ps0, pf0, zSh0, zSf0, zBh0, zBf0, gap, pssi2];

% Define the function Phi, which is the lower block of h(x).
% This block is the expected value of the exogenous state variables.
Phi = [rho_y*ey;rho_y*ey_s;rho_d*d;rho_d*d_s;rho_q*q;rho_q*q_s;rho_eps*eps];

% Define eta as in Schmitt-Grohe and Uribe (2004):
%We have 6 zero-mean shocks
eta=[zeros(4,6);ty 0 0 0 0 0;0 ty 0 0 0 0;0 0 td 0 0 0;0 0 0 td 0 0;0 0 0 0 tq 0;0 0 0 0 0 tq;0 0 0 0 0 0];

%Model equations

%budget constraint HOME
f1 = P0*C0*exp(pi+c) - ( exp(ey)*( 1 - ( exp(d0+d)/(1+exp(d0+d)) )  ) ) ...
- (Sh)*( zSh0*exp(zSh) + exp(ey)*(exp(d0+d)/(1+exp(d0+d)) )) ...
- (Sf)*( zSf0*exp(zSf) + pf0*exp(pf+ey_s)*(exp(d0+d_s)/(1+exp(d0+d_s)) ) ) ... 
- P0*Bh*exp(piu) - Ps0*Bf*exp(piu_s) ... 
+  zSh0*(Shp)*exp(zSh) + zSf0*(Sfp)*exp(zSf)  + zBh0*Bhp*exp(zBh) + zBf0*Bfp*exp(zBf) ; 

%budget constraint FOREIGN
%Note that by imposing market-clearing in financial markets, we have:
%Sh_s = 1 - Sh, Sf_s = 1 - Sf, Bh_s = -Bh, Bf_s = -Bf.
f2 = Ps0*Cs0*exp(pi_s+c_s) - ( pf0*exp(ey_s+pf)*( 1 - ( exp(d0+d_s)/(1+exp(d0+d_s)) )  ) ) ...
- (1 - Sh)*( zSh0*exp(zSh) + exp(ey)*(exp(d0+d)/(1+exp(d0+d)) )) ...
- (1 - Sf)*( zSf0*exp(zSf) + pf0*exp(pf+ey_s)*(exp(d0+d_s)/(1+exp(d0+d_s)) )) ...
- (P0*(-Bh)*exp(piu) + Ps0*(-Bf)*exp(piu_s))  ... 
+ zSh0*(1-Shp)*exp(zSh) + zSf0*(1-Sfp)*exp(zSf) + zBh0*(-Bhp)*exp(zBh) + zBf0*(-Bfp)*exp(zBf);  

%Price index HOME
f3 = (P0^(1-phi))*exp(pi*(1-phi)) ...
    - (    alpa*exp(-(1-phi)*q) + (1-alpa)*(pf0^(1-phi))*exp((pf-q_s)*(1-phi)));
%Price index FOREIGN
f4 = (Ps0^(1-phi))*exp(pi_s*(1-phi)) - ...
    ((1-alpa)*exp(-(1-phi)*q) +     alpa*(pf0^(1-phi))*exp((pf-q_s)*(1-phi)));

%Non-adjusted price index HOME 
f5 = (P0^(1-phi))*exp(piu*(1-phi)) - ...
    (    alpa + (1-alpa)*(pf0^(1-phi))*exp((pf)*(1-phi)));
%Non-adjusted price index FOREIGN
f6 = (Ps0^(1-phi))*exp(piu_s*(1-phi)) - ...
    ((1-alpa) +     alpa*(pf0^(1-phi))*exp((pf)*(1-phi)));

%Market-clearing HOME good;
f7 = alpa*((P0^phi)*C0)*exp(pi*phi + c) ...
    + (1-alpa)*((Ps0^phi)*Cs0)*exp(pi_s*phi + c_s) - exp(ey-(phi-1)*q);

%Euler equation Home equity, HOME
f8 = betta*(exp(-kappa*c))*( zSh0*exp(zShp) + exp(eyp)*( exp(d0+dp)/(1+exp(d0+dp)) ) )*exp( -pip-(gama+gap)*cp ) ...
- zSh0*exp(zSh - pi - (gama+gap)*c)*(1 + (1-eps)*pssi*(Shp-a));
%Euler equation foreign equity, HOME
f9 = betta*(exp(-kappa*c))*( zSf0*exp(zSfp) + pf0*exp(pfp+ey_sp)*(exp(d0+d_sp)/(1+exp(d0+d_sp)) ) )*exp( -pip-(gama+gap)*cp ) ...
- zSf0*exp(zSf - pi - (gama+gap)*c)*(1 + (1-eps)*pssi*(Sfp-af));

%Euler equation Home equity, FOREIGN
f10 = betta*(exp(-kappa*c_s))*( zSh0*exp(zShp) + exp(eyp)*(exp(d0+dp)/(1+exp(d0+dp)) ) )*exp( -pi_sp-(gama-gap)*c_sp ) ...
- zSh0*exp(zSh - pi_s - (gama-gap)*c_s)*(1 - (1-eps)*pssi*(Shp-a));
%Euler equation foreign equity, FOREIGN
f11 = betta*(exp(-kappa*c_s))*( zSf0*exp(zSfp) + pf0*exp(pfp+ey_sp)*(exp(d0+d_sp)/(1+exp(d0+d_sp)) ) )*exp( -pi_sp-(gama-gap)*c_sp ) ...
- zSf0*exp(zSf - pi_s - (gama-gap)*c_s)*(1 - (1-eps)*pssi*(Sfp-af));

%Euler equation Home bond, HOME
f12 = betta*(exp(-kappa*c))*P0*exp(piup  -pip -(gama+gap)*cp) ...
    - zBh0*exp(zBh-pi-(gama+gap)*c)*(1 + (1-eps)*pssi2*(Bhp-b));
%Euler equation foreign bond, HOME
f13 = betta*(exp(-kappa*c))*Ps0*exp(piu_sp-pip -(gama+gap)*cp) ...
    - zBf0*exp(zBf-pi-(gama+gap)*c)*(1 + (1-eps)*pssi2*(Bfp-bf));

%Euler equation Home bond, FOREIGN
f14 = betta*(exp(-kappa*c_s))*P0*exp(piup-  pi_sp-(gama-gap)*c_sp) ...
    - zBh0*exp(zBh-pi_s-(gama-gap)*c_s)*(1 - (1-eps)*pssi2*(Bhp-b));
%Euler equation foreign bond, FOREIGN
f15 = betta*(exp(-kappa*c_s))*Ps0*exp(piu_sp-pi_sp-(gama-gap)*c_sp) ...
    - zBf0*exp(zBf-pi_s-(gama-gap)*c_s)*(1 - (1-eps)*pssi2*(Bfp-bf));



f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15]; % all conditions

% Choose approximation order

approx=3; % third order solution.

% call differentiate_dsge

model=differentiate_dsge(f,yp,y,xp,x,symparams,approx,Phi);
save('model')