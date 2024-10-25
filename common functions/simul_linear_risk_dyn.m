function [yt,xt]=simul_linear_risk_dyn(x0,shocks,mdr,yss,eps_ind,stf,stl,h)

% This function generates simulations using linear risk-corrected policy
% rules computed with Dynare, centered at the approximated SSS. 
%
% It simulates the model from the initial state x0. shocks is a matrix with
% n_e rows and T columns, where n_e is the number of shocks, and T is the
% length of the simulation. mdr is a structure array with Dynare's solution,
% yss the DSS, eps_ind the index of epsilon in the vector of states, stf and
% stl the indices of the first and last state in the vector of Dynare's
% endogenous variables, % (the assumption is that all endogenous variables 
% from y(stf) to y(stl) are states). Finally h is the step size for numerical 
% differentiation.

% The outputs are simulations yt (all the endogenous variables) 
% and xt (all the states) for T+1 periods. 
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[jac_ghx,jac_ghu] = my_jac_yt_dyn(mdr,yss,eps_ind,h);

[aa_ghx,bb_ghx]=size(mdr.ghx);
[~,bb_ghu]=size(mdr.ghu);
lx = bb_ghx;
lu = bb_ghu;

x_dif0 = zeros(lx,1);
x_dif0(eps_ind) = 1;
ysss = dr_yt(mdr,yss,3,x_dif0,zeros(lu,1));
xsss = yss(stf:stl,:);

T=size(shocks,2);
shocks=[zeros(lu,1),shocks,zeros(lu,1)];

xt_f=zeros(lx,T+2);
yt=zeros(aa_ghx,T+2);
xt_f(:,1)=x0-xsss;


for t=1:T+1

    x_f=xt_f(:,t);
    yt(:,t+1) = jac_ghx*x_f + jac_ghu*shocks(:,t+1);
    xt_f(:,t+1) = yt(stf:stl,t+1);

end

yt=yt(:,1:T+1);
xt=xt_f(:,1:T+1);


yt=yt+repmat(ysss,1,T+1);
xt=xt(:,:)+repmat(xsss,1,T+1);