function y = system_soe_o2(x,sol1,P)
%SOE with PAC as auxiliary model
%Full version with income and interest rate shocks
%One-parameter perturbation
%This function evaluates second-order derivatives of Euler w.r.t. states

%parameters
betta = P(1);
gama = P(2);
rho_y = P(3);
rho_r = P(4);
ty = P(5);
tr = P(6);
R0 = P(7);
b_bar = P(8);
C0 = P(9);
psi1 = P(10);
psi2 = P(11);
uy_ur_corr = P(12);

%First-order solution
hb = sol1(1);
hr = sol1(2);
hy = sol1(3);

%second-order coeffs
hs2 = x(1);
hb2 = x(2);
hr2 = x(3);
hy2 = x(4);
hbr = x(5);
hby = x(6);
hry = x(7);

%Equations
Ub2 = betta*gama*((hb^2*(psi1*hb^2 + hb2/R0))/C0 - (hb2*(R0 - hb))/(C0*R0)) - hb2*psi1 - (gama*(psi1*hb^2 + hb2/R0))/(C0*R0) - (gama*(R0 - hb)^2*(gama + 1))/(C0^2*R0^3) + (2*gama*hb*psi1*(R0 - hb))/(C0*R0) + (betta*gama*hb^2*(R0 - hb)^2*(gama + 1))/(C0^2*R0^2);
Ur2 = betta*gama*((rho_r^2*(psi1*hr^2 - (2*hr)/R0 + b_bar/R0 + hr2/R0))/C0 + (hr^2*(psi1*hb^2 + hb2/R0))/C0 - (hr2*(R0 - hb))/(C0*R0) + (2*hr*rho_r*(hbr/R0 - hb/R0 + hb*hr*psi1))/C0) - 1/R0 - hr2*psi1 - (gama*(psi1*hr^2 - (2*hr)/R0 + b_bar/R0 + hr2/R0))/(C0*R0) + betta*gama*((hr*(R0 - hb))/(C0*R0) + (rho_r*(b_bar - hr))/(C0*R0))^2*(gama + 1) + (2*gama*(b_bar - hr)*(hr*psi1 - 1/R0))/(C0*R0) - (gama*(b_bar - hr)^2*(gama + 1))/(C0^2*R0^3);
Uy2 = betta*gama*((hy^2*(psi1*hb^2 + hb2/R0))/C0 + (rho_y^2*(psi1*hy^2 + hy2/R0 - 1))/C0 + (2*hy*rho_y*(hby/R0 + hb*hy*psi1))/C0 - (hy2*(R0 - hb))/(C0*R0)) - hy2*psi1 - (gama*(psi1*hy^2 + hy2/R0 - 1))/(C0*R0) + betta*gama*((hy*(R0 - hb))/(C0*R0) + (rho_y*(R0 - hy))/(C0*R0))^2*(gama + 1) - (gama*(R0 - hy)^2*(gama + 1))/(C0^2*R0^3) + (2*gama*hy*psi1*(R0 - hy))/(C0*R0);
Ubr = betta*gama*((hb*hr*(psi1*hb^2 + hb2/R0))/C0 - (hbr*(R0 - hb))/(C0*R0) + (hb*rho_r*(hbr/R0 - hb/R0 + hb*hr*psi1))/C0) - hbr*psi1 - (gama*(hbr/R0 - hb/R0 + hb*hr*psi1))/(C0*R0) + (gama*(R0 - hb)*(hr*psi1 - 1/R0))/(C0*R0) + (gama*hb*psi1*(b_bar - hr))/(C0*R0) - (gama*(R0 - hb)*(b_bar - hr)*(gama + 1))/(C0^2*R0^3) + (betta*gama*hb*(R0 - hb)*((hr*(R0 - hb))/(C0*R0) + (rho_r*(b_bar - hr))/(C0*R0))*(gama + 1))/(C0*R0);
Uby = betta*gama*((hb*hy*(psi1*hb^2 + hb2/R0))/C0 + (hb*rho_y*(hby/R0 + hb*hy*psi1))/C0 - (hby*(R0 - hb))/(C0*R0)) - hby*psi1 - (gama*(hby/R0 + hb*hy*psi1))/(C0*R0) + (gama*hb*psi1*(R0 - hy))/(C0*R0) + (gama*hy*psi1*(R0 - hb))/(C0*R0) - (gama*(R0 - hb)*(R0 - hy)*(gama + 1))/(C0^2*R0^3) + (betta*gama*hb*(R0 - hb)*((hy*(R0 - hb))/(C0*R0) + (rho_y*(R0 - hy))/(C0*R0))*(gama + 1))/(C0*R0);
Ury = betta*gama*((hr*hy*(psi1*hb^2 + hb2/R0))/C0 + (hr*rho_y*(hby/R0 + hb*hy*psi1))/C0 - (hry*(R0 - hb))/(C0*R0) + (hy*rho_r*(hbr/R0 - hb/R0 + hb*hr*psi1))/C0 + (rho_r*rho_y*(hry/R0 - hy/R0 + hr*hy*psi1))/C0) - hry*psi1 - (gama*(hry/R0 - hy/R0 + hr*hy*psi1))/(C0*R0) + betta*gama*((hy*(R0 - hb))/(C0*R0) + (rho_y*(R0 - hy))/(C0*R0))*((hr*(R0 - hb))/(C0*R0) + (rho_r*(b_bar - hr))/(C0*R0))*(gama + 1) + (gama*(R0 - hy)*(hr*psi1 - 1/R0))/(C0*R0) + (gama*hy*psi1*(b_bar - hr))/(C0*R0) - (gama*(R0 - hy)*(b_bar - hr)*(gama + 1))/(C0^2*R0^3);
Us2 = uy_ur_corr*((2*betta*gama*tr*ty*(hry/R0 - hy/R0 + hr*hy*psi1))/C0 + (2*betta*gama*tr*ty*(R0 - hy)*(b_bar - hr)*(gama + 1))/(C0^2*R0^2)) + 2*betta*psi2 - hs2*psi1 + betta*gama*(hs2/(C0*R0) - (hs2*(R0 - hb))/(C0*R0)) - (gama*hs2)/(C0*R0^2) + (betta*gama*ty^2*(psi1*hy^2 + hy2/R0 - 1))/C0 + (betta*gama*tr^2*(psi1*hr^2 - (2*hr)/R0 + b_bar/R0 + hr2/R0))/C0 + (betta*gama*tr^2*(b_bar - hr)^2*(gama + 1))/(C0^2*R0^2) + (betta*gama*ty^2*(R0 - hy)^2*(gama + 1))/(C0^2*R0^2);
 
y = [Ub2;Ur2;Uy2;Ubr;Uby;Ury;Us2];