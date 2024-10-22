function y = sol_soe_o1_only_y(P)
%SOE with PAC as auxiliary model
%Basic version with only income shocks
%One-parameter perturbation
%This function computes the closed-form first-order solution

betta = P(1);
gama = P(2);
rho_y = P(3);
ty = P(5);
R0 = P(7);
b_bar = P(8);
C0 = P(9);
psi1 = P(10);



hb =(gama - (C0^2*R0^4*psi1^2 + 2*C0*R0^4*betta*gama*psi1 + 2*C0*R0^2*gama*psi1 + R0^4*betta^2*gama^2 - 2*R0^2*betta*gama^2 + gama^2)^(1/2) + C0*R0^2*psi1 + R0^2*betta*gama)/(2*R0*betta*gama);
hy =(gama/(C0*R0) - (betta*gama*rho_y)/C0)/(psi1 + gama/(C0*R0^2) - betta*gama*(rho_y/(C0*R0) - (R0 - hb)/(C0*R0)));
 
y= [hb hy];
 