%SOE with PAC as auxiliary model
%Full version with income and interest rate shocks
%One-parameter perturbation
%This script computes analytical derivatives up to third order

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

syms hb3 hr3 hy3;
syms gb3 gr3 gy3;
syms hb2r hbr2 hb2y hby2 hr2y hry2;
syms gb2r gbr2 gb2y gby2 gr2y gry2;
syms hbry gbry;
syms hs2b hs2r hs2y;
syms gs2b gs2r gs2y;

%Tomorrow's bonds
btp = b_bar + hb*(b_1-b_bar) + hr*rt + hy*yt ...
    + 0.5*( hb2*(b_1-b_bar)^2 + hr2*rt^2 + hy2*yt^2 + hs2*s^2) ...
    + hbr*(b_1-b_bar)*rt + hby*(b_1-b_bar)*yt+ hry*rt*yt ...
    + (1/6)*( hb3*(b_1-b_bar)^3 + hr3*rt^3 + hy3*yt^3 ) ...
    + 0.5*( hb2r*rt*(b_1-b_bar)^2 + hbr2*(b_1-b_bar)*rt^2 + hb2y*yt*(b_1-b_bar)^2 + hby2*(b_1-b_bar)*yt^2 + hr2y*yt*rt^2 + hry2*rt*yt^2 + hs2b*(b_1-b_bar)*s^2 + hs2r*rt*s^2 + hs2y*yt*s^2 ) ...
    + hbry*(b_1-b_bar)*rt*yt;

%Today's consumption
ct = 1 + gb*(b_1-b_bar) + gr*rt + gy*yt ...
    + 0.5*( gb2*(b_1-b_bar)^2 + gr2*rt^2 + gy2*yt^2 + gs2*s^2) ...
    + gbr*(b_1-b_bar)*rt + gby*(b_1-b_bar)*yt+ gry*rt*yt ...
    + (1/6)*( gb3*(b_1-b_bar)^3 + gr3*rt^3 + gy3*yt^3 ) ...
    + 0.5*( gb2r*rt*(b_1-b_bar)^2 + gbr2*(b_1-b_bar)*rt^2 + gb2y*yt*(b_1-b_bar)^2 + gby2*(b_1-b_bar)*yt^2 + gr2y*yt*rt^2 + gry2*rt*yt^2 + gs2b*(b_1-b_bar)*s^2 + gs2r*rt*s^2 + gs2y*yt*s^2 ) ...
    + gbry*(b_1-b_bar)*rt*yt;

%Tomorrow's exo states (log output and log interest rate)
yt1 = rho_y*yt + s*ty*uy;
rt1 = rho_r*rt + s*tr*ur;

%Tomorrow's consumption
ct1 = 1 + gb*(btp-b_bar) + gr*rt1 + gy*yt1 ...
    + 0.5*( gb2*(btp-b_bar)^2 + gr2*rt1^2 + gy2*yt1^2 + gs2*s^2) ...
    + gbr*(btp-b_bar)*rt1 + gby*(btp-b_bar)*yt1+ gry*rt1*yt1 ...
    + (1/6)*( gb3*(btp-b_bar)^3 + gr3*rt1^3 + gy3*yt1^3 ) ...
    + 0.5*( gb2r*rt1*(btp-b_bar)^2 + gbr2*(btp-b_bar)*rt1^2 + gb2y*yt1*(btp-b_bar)^2 + gby2*(btp-b_bar)*yt1^2 + gr2y*yt1*rt1^2 + gry2*rt1*yt1^2 + gs2b*(btp-b_bar)*s^2 + gs2r*rt1*s^2 + gs2y*yt1*s^2 ) ...
    + gbry*(btp-b_bar)*rt1*yt1;

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

%New order 3 derivatives for Euler
%gb3 gr3 gy3
Ub3 = diff(Ub2,b_1);
Ur3 = diff(Ur2,rt);
Uy3 = diff(Uy2,yt);

%gb2r gbr2 gb2y gby2 gr2y gry2
Ub2r = diff(Ubr,b_1);
Ubr2 = diff(Ubr,rt);
Ub2y = diff(Uby,b_1);
Uby2 = diff(Uby,yt);
Ur2y = diff(Ury,rt);
Ury2 = diff(Ury,yt);

%hs2b hs2r hs2y
Us2b = diff(Us2,b_1);
Us2r = diff(Us2,rt);
Us2y = diff(Us2,yt);

%hbry
Ubry = diff(Ubr,yt);

%New order 3 derivatives for BC
BCb3 = diff(BCb2,b_1);
BCr3 = diff(BCr2,rt);
BCy3 = diff(BCy2,yt);

%gb2r gbr2 gb2y gby2 gr2y gry2

BCb2r = diff(BCbr,b_1);
BCbr2 = diff(BCbr,rt);
BCb2y = diff(BCby,b_1);
BCby2 = diff(BCby,yt);
BCr2y = diff(BCry,rt);
BCry2 = diff(BCry,yt);

%hs2b hs2r hs2y
BCs2b = diff(BCs2,b_1);
BCs2r = diff(BCs2,rt);
BCs2y = diff(BCs2,yt);

BCbry = diff(BCbr,yt);

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

BCb3 = subs(BCb3);
gb3 = solve(BCb3,gb3)
BCr3 = subs(BCr3);
gr3 = solve(BCr3,gr3)
BCy3 = subs(BCy3);
gy3 = solve(BCy3,gy3)

%gb2r gbr2 gb2y gby2 gr2y gry2

BCb2r = subs(BCb2r);
gb2r = solve(BCb2r,gb2r)
BCbr2 = subs(BCbr2);
gbr2 = solve(BCbr2,gbr2)
BCb2y = subs(BCb2y);
gb2y = solve(BCb2y,gb2y)
BCby2 = subs(BCby2);
gby2 = solve(BCby2,gby2)
BCr2y = subs(BCr2y);
gr2y = solve(BCr2y,gr2y)
BCry2 = subs(BCry2);
gry2 = solve(BCry2,gry2)

%hs2b hs2r hs2y
BCs2b = subs(BCs2b);
gs2b = solve(BCs2b,gs2b)
BCs2r = subs(BCs2r);
gs2r = solve(BCs2r,gs2r)
BCs2y = subs(BCs2y);
gs2y = solve(BCs2y,gs2y)

BCbry = subs(BCbry);
gbry = solve(BCbry,gbry)
%------------------------------------------
%Substitute in Euler derivatives (only order 3)  

Ub3 = subs(Ub3)
Ur3 = subs(Ur3)
Uy3 = subs(Uy3)

Ub2r = subs(Ub2r)
Ubr2 = subs(Ubr2)
Ub2y = subs(Ub2y)
Uby2 = subs(Uby2)
Ur2y = subs(Ur2y)
Ury2 = subs(Ury2)

Ubry = subs(Ubry)

%Special treatment for Us2b, Us2r and Us2y (depend on innovations)
Us2b = subs(Us2b);
Us2b = collect(Us2b,[ur uy]);
Us2b = subs(Us2b,ur^2,1); %ur^2 for var(ur)
Us2b = subs(Us2b,uy^2,1); %uy^2 for var(uy)
Us2b = subs(Us2b,ur*uy,uy_ur_corr); %ur*uy for cov(ur,uy)

Us2b

Us2r = subs(Us2r);
Us2r = collect(Us2r,[ur uy]);
Us2r = subs(Us2r,ur^2,1);
Us2r = subs(Us2r,uy^2,1);
Us2r = subs(Us2r,ur*uy,uy_ur_corr);

Us2r

Us2y = subs(Us2y);
Us2y = collect(Us2y,[ur uy]);
Us2y = subs(Us2y,ur^2,1);
Us2y = subs(Us2y,uy^2,1);
Us2y = subs(Us2y,ur*uy,uy_ur_corr);

Us2y


