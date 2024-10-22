function y = sol_soe_o2_only_y(sol1,P)
%SOE with PAC as auxiliary model
%Basic version with only income shocks
%One-parameter perturbation
%This function computes the closed-form second-order solution

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

hb = sol1(1);
hy = sol1(2);

hb2 =((betta*gama*hb^4*psi1)/C0 - (gama*hb^2*psi1)/(C0*R0) - (gama*(R0 - hb)^2*(gama + 1))/(C0^2*R0^3) + (2*gama*hb*psi1*(R0 - hb))/(C0*R0) + (betta*gama*hb^2*(R0 - hb)^2*(gama + 1))/(C0^2*R0^2))/(psi1 + gama/(C0*R0^2) + betta*gama*((R0 - hb)/(C0*R0) - hb^2/(C0*R0)));
hby =(betta*gama*((hb*hy*(psi1*hb^2 + hb2/R0))/C0 + (hb^2*hy*psi1*rho_y)/C0) - (gama*hb*hy*psi1)/(C0*R0) + (gama*hb*psi1*(R0 - hy))/(C0*R0) + (gama*hy*psi1*(R0 - hb))/(C0*R0) - (gama*(R0 - hb)*(R0 - hy)*(gama + 1))/(C0^2*R0^3) + (betta*gama*hb*(R0 - hb)*((hy*(R0 - hb))/(C0*R0) + (rho_y*(R0 - hy))/(C0*R0))*(gama + 1))/(C0*R0))/(psi1 + gama/(C0*R0^2) + betta*gama*((R0 - hb)/(C0*R0) - (hb*rho_y)/(C0*R0)));
hy2 =(betta*gama*((hy^2*(psi1*hb^2 + hb2/R0))/C0 + (rho_y^2*(psi1*hy^2 - 1))/C0 + (2*hy*rho_y*(hby/R0 + hb*hy*psi1))/C0) - (gama*(psi1*hy^2 - 1))/(C0*R0) + betta*gama*((hy*(R0 - hb))/(C0*R0) + (rho_y*(R0 - hy))/(C0*R0))^2*(gama + 1) - (gama*(R0 - hy)^2*(gama + 1))/(C0^2*R0^3) + (2*gama*hy*psi1*(R0 - hy))/(C0*R0))/(psi1 + gama/(C0*R0^2) + betta*gama*(- rho_y^2/(C0*R0) + (R0 - hb)/(C0*R0)));
hs2 =(2*betta*psi2 + (betta*gama*ty^2*(psi1*hy^2 + hy2/R0 - 1))/C0 + (betta*gama*ty^2*(R0 - hy)^2*(gama + 1))/(C0^2*R0^2))/(psi1 + gama/(C0*R0^2) - betta*gama*(1/(C0*R0) - (R0 - hb)/(C0*R0)));


y = [hs2 hb2 hb2 hy2];