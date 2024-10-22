%SOE with PAC as auxiliary model
%Full version with income and interest rate shocks
%One-parameter perturbation
%This script computes analytical derivatives up to second order

clear;clc;

%parameters
syms betta gama rho_y rho_r rho_eps ty tr R0 b_bar C0 psi1 psi2 uy_ur_corr real
%states (s is the parturbation parameter)
syms b_1 yt rt s real
%innovations
syms uy ur;
%Coefficients of decision rules, h for bonds, g for consumption
syms hb hr hy;
syms gb gr gy;
syms hb2 hr2 hy2 hbr hby hry;
syms gb2 gr2 gy2 gbr gby gry;
syms gs2 hs2;

%Today's consumption
ct = 1 + gb*(b_1-b_bar) + gr*rt + gy*yt + 0.5*( gb2*(b_1-b_bar)^2 + gr2*rt^2 + gy2*yt^2 + gs2*s^2) + gbr*(b_1-b_bar)*rt + gby*(b_1-b_bar)*yt+ gry*rt*yt;
%Tomorrow's bonds
btp = b_bar + hb*(b_1-b_bar) + hr*rt + hy*yt + 0.5*( hb2*(b_1-b_bar)^2 + hr2*rt^2 + hy2*yt^2 + hs2*s^2) + hbr*(b_1-b_bar)*rt + hby*(b_1-b_bar)*yt+ hry*rt*yt;
%Tomorrow's exo states (log output and log interest rate)
yt1 = rho_y*yt + s*ty*uy;
rt1 = rho_r*rt + s*tr*ur;
%Tomorrow's consumption
ct1 = 1 + gb*(btp-b_bar) + gr*rt1 + gy*yt1 + 0.5*( gb2*(btp-b_bar)^2 + gr2*rt1^2 + gy2*yt1^2 + gs2*s^2) + gbr*(btp-b_bar)*rt1 + gby*(btp-b_bar)*yt1+ gry*rt1*yt1;

%budget constraint
BC =   exp(yt) + b_1 ...
    - C0*(ct) - ((R0*exp(rt))^(-1))*btp - 0.5*psi1*(1-s^2)*(btp-b_bar)^2;
%Euler equation
U = (betta*(1+psi2*(s^2)))*(ct1^-gama) ...
    - (ct^-gama)*( 1/(R0*exp(rt)) + psi1*(1-s^2)*(btp-b_bar) );

%Derivatives of Euler w.r.t states
Ub = diff(U,b_1);
Ur = diff(U,rt);
Uy = diff(U,yt);
Ub2 = diff(Ub,b_1);
Ur2 = diff(Ur,rt);
Uy2 = diff(Uy,yt);
Ubr = diff(Ub,rt);
Uby = diff(Ub,yt);
Ury = diff(Ur,yt);
Us2 = diff(diff(U,s),s);

%Derivatives of BC w.r.t states
BCb = diff(BC,b_1);
BCr = diff(BC,rt);
BCy = diff(BC,yt);
BCb2 = diff(BCb,b_1);
BCr2 = diff(BCr,rt);
BCy2 = diff(BCy,yt);
BCbr = diff(BCb,rt);
BCby = diff(BCb,yt);
BCry = diff(BCr,yt);
BCs2 = diff(diff(BC,s),s);

%Evaluate at DSS
b_1 = b_bar; rt = 0; yt = 0; s = 0; 

%Substitute in BC derivatives and solve for g coeffs and functions of h coeffs and params
BCb = subs(BCb);
gb=solve(BCb,gb)
BCr = subs(BCr);
gr=solve(BCr,gr)
BCy = subs(BCy);
gy=solve(BCy,gy)
BCb2 = subs(BCb2);
gb2=solve(BCb2,gb2)
BCr2 = subs(BCr2);
gr2=solve(BCr2,gr2)
BCy2 = subs(BCy2);
gy2=solve(BCy2,gy2)
BCbr = subs(BCbr);
gbr=solve(BCbr,gbr)
BCby = subs(BCby);
gby=solve(BCby,gby)
BCry = subs(BCry);
gry=solve(BCry,gry)
BCs2 = subs(BCs2);
gs2=solve(BCs2,gs2)

%Substitute in Euler derivatives 
U = subs(U)
Ub = subs(Ub)
Ur = subs(Ur)
Uy = subs(Uy)
Ub2 = subs(Ub2)
Ur2 = subs(Ur2)
Uy2 = subs(Uy2)
Ubr = subs(Ubr)
Uby = subs(Uby)
Ury = subs(Ury)

%Special treatment for Us2 (depends on innovations)
Us2 = subs(Us2);
Us2 = collect(Us2,[ur uy]);

Us2 = subs(Us2,ur^2,1); %ur^2 for var(ur)
Us2 = subs(Us2,uy^2,1); %uy^2 for var(uy)
Us2 = subs(Us2,ur*uy,uy_ur_corr);  %ur*uy for cov(ur,uy)

Us2

