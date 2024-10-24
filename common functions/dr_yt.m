function y = dr_yt(mdr,yss,approx,x,u)
%This function evaluates decision rules of a DSGE model solved with Dynare.
%The order of approximation can be either 2 or 3.
%Inputs are:
%mdr: struct with derivatives of decision rules, generated by Dynare.
%yss: (nx+ny)-by-1 vector deterministic steady-state
%approx: order of approximation
%x: nx-by-1 vector of state variables
%u: nu-by-N matrix of innovations, where N>=1, and each column is a particular
%realization of innovations.
%
%The output is a (nx+ny)-by-N matrix, where each column represents
%time t model variables under a particular realization of innovations.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

ghx = mdr.ghx;
ghu = mdr.ghu;
ghxx = mdr.ghxx;
ghxu = mdr.ghxu;
ghuu = mdr.ghuu;
ghs2 = mdr.ghs2;
if approx>2
    ghxxx = mdr.ghxxx;
    ghxxu = mdr.ghxxu;
    ghxuu = mdr.ghxuu;
    ghuuu = mdr.ghuuu;
    ghxss = mdr.ghxss;
    ghuss = mdr.ghuss;
end

kx2 = mykron1(x,x);
ku2 = mykron1(u,u);
kxu = mykron1(x,u);
kx3 = mykron1(kx2,x);
kx2u = mykron1(kx2,u);
kxu2 = mykron1(x,ku2);
ku3 = mykron1(ku2,u);


yp = yss + 0.5*ghs2 + ghx*x + ghu*u + 0.5*ghxx*kx2 + 0.5*ghuu*ku2 + ghxu*kxu;

if approx>2
    
    yp = yp + (1/6)*ghxxx*kx3 ...
        + 0.5*ghxxu*kx2u ...
        + 0.5*ghxuu*kxu2 ...
        + (1/6)*ghuuu*ku3 ...
        + 0.5*ghxss*x + 0.5*ghuss*u;
    
end

y = yp;

