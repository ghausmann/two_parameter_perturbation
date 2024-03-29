function g = dr_gt(mdr,gss,approx,x)
%This function evaluates decision rules for control variables of a DSGE
%model solved with the algorithm by Levintal (2017).
%The order of approximation can go from 2 up to 5.
%Inputs are:
%mdr: struct with derivatives of decision rules, generated by the algorithm.
%gss: ny-by-1 vector of deterministic steady-state for controls.
%approx: order of approximation
%x: nx-by-N vector of state variables, where N>=1, and each column is a particular
%realization of states.
%
%The output is a ny-by-N matrix, where each column represents
%time t control variables under a particular realization of states.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

[~,s2] = size(x);

x = [x;ones(1,s2)];

gx = mdr.gx;
gx2 = mdr.gxx;
kx2 = mykron1(x,x);
if approx>2
    gx3 = mdr.gxxx;
    kx3 = mykron1(kx2,x);
end
if approx>3
    gx4 = mdr.gxxxx;
    kx4 = mykron1(kx3,x);
end
if approx>4
    gx5 = mdr.gxxxxx;
    kx5 = mykron1(kx4,x);
end

if approx==2
    g = gss +  gx*x + 0.5*gx2*kx2;
end
if approx==3
    g = gss +  gx*x + 0.5*gx2*kx2 + (1/6)*gx3*kx3;
end
if approx==4
    g = gss +  gx*x + 0.5*gx2*kx2 + (1/6)*gx3*kx3 + (1/24)*gx4*kx4;
end
if approx==5
    g = gss +  gx*x + 0.5*gx2*kx2 + (1/6)*gx3*kx3 + (1/24)*gx4*kx4 + (1/120)*gx5*kx5;
end

