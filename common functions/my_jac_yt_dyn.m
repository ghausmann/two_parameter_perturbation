function [jac_ghx,jac_ghu] = my_jac_yt_dyn(mdr,yss,eps_ind,h)
%This function computes the jacobian of the policy rules for all Dynare's
%endoganous variables, fixing sigma=epsilon=1.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[aa_ghx,bb_ghx]=size(mdr.ghx);
[aa_ghu,bb_ghu]=size(mdr.ghu);

lx = bb_ghx;
lu = bb_ghu;

ghx_jac = zeros(aa_ghx,bb_ghx);
ghu_jac = zeros(aa_ghu,bb_ghu);

x_dif0 = zeros(lx,1);
x_dif0(eps_ind) = 1;
ysss = dr_yt(mdr,yss,3,x_dif0,zeros(lu,1));

for n = 1:lx

    x_dif1 = x_dif0;
    x_dif1(n) = x_dif1(n) + h;
    y1 = dr_yt(mdr,yss,3,x_dif1,zeros(lu,1));
    ghx_jac(:,n) = (y1-ysss)/h;

end

ghx_jac(:,eps_ind) = zeros(aa_ghx,1);

for n = 1:lu

    u_dif = zeros(lu,1);
    u_dif(n) = h;
    y1 = dr_yt(mdr,yss,3,x_dif0,u_dif);
    ghu_jac(:,n) = (y1-ysss)/h;

end

jac_ghx = ghx_jac;
jac_ghu = ghu_jac;