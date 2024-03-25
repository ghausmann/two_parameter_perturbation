function [yt,xt]=simul_mod_pruning3(x0,shocks,nyss,nxss,eta,derivs)
% [yt,xt]=simul_mod_pruning3(x0,shocks,nyss,nxss,eta,derivs) is a
% heavily modified version of the original function simul.m by Oren
% Levintal (2016). You can use it to simulate a two-parameter perturbation
% model where third-order policy functions are approximated around the
% stochastic steady-state of the state variables.
%
%It simulates the model from the initial state x0. shocks is a matrix with
%n_e rows and T columns, where n_e is the number of shocks (corresponds to
%the columns of eta), and T is the length of the simulation. The function
% returns yt and xt for T+1 periods. The first period is the initial state
% and the rest T periods correspond to the shocks. The pruning algorithm is
% a version of Andreasen, Fernandez-Villaverde and Rubio-Ramirez (2018)
% "The Pruned State-Space System for Non-Linear DSGE Models: Theory and
% Empirical Applications".
%
% KEY MODIFICATION: This code treats the perturbation objects (sigma and epsilon) as
% parameters). To do so:
% a) It computes the jacobians of third-order policy functions
% (fixing sigma=1 and epsilon=1) with respect to the states.
% b) This new jacobians incorporate up to third-order perturbation
% effects, and replace the jacobians of the auxiliary deterministic model (the
% first-order solution).
% c) It takes a single evaluation of the policy functions (fixing sigma=1
% and epsilon=1) to compute the stochastic steady-state of the controls,
% which replaces their DSS values in the pruned functions.
% d) It then constructs the pruned functions as in Andreasen,
% Fernandez-Villaverde and Rubio-Ramirez (2018), imposing that perturbation
% variables equal 0 (sigma_t=0, epsilon_t=0) to avoid double accounting of
% perturbation effects.
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

% Compute jacobians accounting for third-order perturbation effects
derivs_hx_c = my_jac_hx(derivs,nxss,3,0.0001);
derivs_gx_c = my_jac_gx(derivs,nxss,nyss,3,0.0001);

% Evaluate policy rules at (sigma,epsilon)=(1,1) to find SSS
nnyss = dr_gt(derivs,nyss,3,(x_0-nxss));
nnxss = dr_ht(derivs,nxss,3,(x_0-nxss));

T=size(shocks,2);
n_y=size(derivs.gx,1);
n_x=size(derivs.hx,1);
n_e=size(shocks,1);
shocks=[zeros(n_e,1),shocks,zeros(n_e,1)];

xt_f=zeros(n_x+1,T+2);
yt=zeros(n_y,T+2);
xt_s=zeros(n_x+1,T+2);
xt_rd=zeros(n_x+1,T+2);
xt_f(1:end-1,1)=x0-nxss;


for t=1:T+1
    x_f=xt_f(:,t);
    xt_f(1:end-1,t+1)=derivs_hx_c*x_f+eta*shocks(:,t+1);
    
    x_s=xt_s(:,t);
    x_f2=kron(x_f,x_f);
    xt_s(1:end-1,t+1)=derivs_hx_c*x_s+derivs.hxx*x_f2/2;
    
    x_rd=xt_rd(:,t);
    x_f3=kron(x_f2,x_f);
    x_f_x_s=kron(x_f,x_s);
    xt_rd(1:end-1,t+1)=derivs_hx_c*x_rd+derivs.hxx*(2*x_f_x_s)/2+derivs.hxxx*x_f3/6;
    
    
    yt(:,t)=derivs_gx_c*(x_f+x_s+x_rd)+derivs.gxx*(x_f2+2*x_f_x_s)/2 ...
        +derivs.gxxx*(x_f3)/6;
    
end

yt=yt(:,1:T+1);
xt=xt_f(:,1:T+1) + xt_s(:,1:T+1) + xt_rd(:,1:T+1);


yt=yt+repmat(nnyss,1,T+1);
xt=xt(1:end-1,:)+repmat(nnxss,1,T+1);


