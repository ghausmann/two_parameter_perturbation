%SOE with PAC as auxiliary model
%Basic version with only income shocks
%One-parameter perturbation
%This script computes analytical derivatives up to second order

clear;clc;

%parameters
syms betta gama rho_y ty R0 b_bar C0 psi1 psi2 real
%states (s is the parturbation parameter)
syms b_1 yt s real
%innovation
syms uy;
%Coefficients of decision rules, h for bonds, g for consumption
syms hb hy;
syms gb gy;
syms hb2 hy2 hby;
syms gb2 gy2 gby;
syms gs2 hs2;


%Today's consumption
ct = 1 + gb*(b_1-b_bar) + gy*yt + 0.5*( gb2*(b_1-b_bar)^2 + gy2*yt^2 + gs2*s^2) + gby*(b_1-b_bar)*yt;
%Tomorrow's bonds
btp = b_bar + hb*(b_1-b_bar) + hy*yt + 0.5*( hb2*(b_1-b_bar)^2 + hy2*yt^2 + hs2*s^2) + hby*(b_1-b_bar)*yt;
%Tomorrow's exo state (log output)
yt1 = rho_y*yt + s*ty*uy;
%Tomorrow's consumption
ct1 = 1 + gb*(btp-b_bar) + gy*yt1 + 0.5*( gb2*(btp-b_bar)^2 + gy2*yt1^2 + gs2*s^2) + gby*(btp-b_bar)*yt1;

%budget constraint
BC =   exp(yt) + b_1 ...
    - C0*(ct) - ((R0)^(-1))*btp - 0.5*psi1*(1-s^2)*(btp-b_bar)^2;
%Euler equation
U = (betta*(1+psi2*(s^2)))*(ct1^-gama) ...
    - (ct^-gama)*( 1/(R0) + psi1*(1-s^2)*(btp-b_bar) );

%Derivatives of Euler w.r.t states
Ub = diff(U,b_1);
Uy = diff(U,yt);
Ub2 = diff(Ub,b_1);
Uy2 = diff(Uy,yt);
Uby = diff(Ub,yt);
Us2 = diff(diff(U,s),s);

%Derivatives of BC w.r.t states
BCb = diff(BC,b_1);
BCy = diff(BC,yt);
BCb2 = diff(BCb,b_1);
BCy2 = diff(BCy,yt);
BCby = diff(BCb,yt);
BCs2 = diff(diff(BC,s),s);

b_1 = b_bar; yt = 0; s = 0; 

%Substitute in BC derivatives and solve for g coeffs and functions of h coeffs and params
BCb = subs(BCb);
gb=solve(BCb,gb)
BCy = subs(BCy);
gy=solve(BCy,gy)
BCb2 = subs(BCb2);
gb2=solve(BCb2,gb2)
BCy2 = subs(BCy2);
gy2=solve(BCy2,gy2)
BCby = subs(BCby);
gby=solve(BCby,gby)
BCs2 = subs(BCs2);
gs2=solve(BCs2,gs2)

%Substitute in Euler derivatives
Ub = subs(Ub)
Uy = subs(Uy)
Ub2 = subs(Ub2)
Uy2 = subs(Uy2)
Uby = subs(Uby)

%Special treatment for Us2 (depends on innovations)
Us2 = subs(Us2);
Us2 = collect(Us2,uy); 
Us2 = subs(Us2,uy^2,1); %uy^2 for var(uy)

%Closed-form solution for 1st order (pick the minus root)
hhb = solve(Ub,hb)
hhy = solve(Uy,hy)
%Closed-form solution for 2nd order
hhb2 = solve(Ub2,hb2)
hhby = solve(Uby,hby)
hhy2 = solve(Uy2,hy2)
hhs2 = solve(Us2,hs2)
