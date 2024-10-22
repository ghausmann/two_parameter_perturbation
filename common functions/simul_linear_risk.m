function [yt,xt]=simul_linear_risk(x0,shocks,nyss,nxss,eta,derivs,approx,h)
% This is a heavily modified version of simul.m by Oren Levintal (2016).
% You can use it to generate simulations from linear risk-corrected policy
% rules obtained by two-parameter perturbation, centered at the
% approximated SSS.
%
%It simulates the model from the initial state x0. shocks is a matrix with
%n_e rows and T columns, where n_e is the number of shocks (corresponds to
%the columns of eta), and T is the length of the simulation. The function
% returns yt and xt for T+1 periods. The first period is the initial state
% and the rest T periods correspond to the shocks. 
%
% CAREFUL! Make sure that the exogenous perturbation variable epsilon is ordered as
% the last state variable.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

lx = numel(nxss);
x0(lx) = 0;
x_0 = nxss;
x_0(lx) = 1;
derivs.hx(lx,lx) = 1;

% Compute jacobians of policy rules after evaluating at (sigma,epsilon)=(1,1) 
derivs_hx_c = my_jac_hx(derivs,nxss,approx,h);
derivs_gx_c = my_jac_gx(derivs,nxss,nyss,approx,h);

% Evaluate policy rules at (sigma,epsilon)=(1,1) to find SSS
nnyss = dr_gt(derivs,nyss,approx,(x_0-nxss));
nnxss = dr_ht(derivs,nxss,approx,(x_0-nxss));

T=size(shocks,2);
n_y=size(derivs.gx,1);
n_x=size(derivs.hx,1);
n_e=size(shocks,1);
shocks=[zeros(n_e,1),shocks,zeros(n_e,1)];

xt_f=zeros(n_x+1,T+2);
yt=zeros(n_y,T+2);
xt_f(1:end-1,1)=x0-nxss;


for t=1:T+1
    x_f=xt_f(:,t);
    xt_f(1:end-1,t+1)=derivs_hx_c*x_f+eta*shocks(:,t+1);
    
    yt(:,t)=derivs_gx_c*(x_f);
    
end

yt=yt(:,1:T+1);
xt=xt_f(:,1:T+1);


yt=yt+repmat(nnyss,1,T+1);
xt=xt(1:end-1,:)+repmat(nnxss,1,T+1);


