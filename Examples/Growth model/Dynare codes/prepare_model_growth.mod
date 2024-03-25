/*
%Neoclassical growth model.
%Auxiliary model (sigma,eps)=(0,0) uses a discount factor different than 
%the one intended for the model of interest (sigma,eps)=(1,1) to obtain
%approximations around a point different than the DSS.
*/

% Variables are k (log capital), z (log TFP), eps (perturbation parameter), 
% c (log consumption), ei (log investment), and ey (log output). 
%Make sure that the declaration order coincides with the dr order, to avoid a mess!
var ei,k,z,eps,c,ey;
varexo u;

parameters betta, gama, alpa, delta, rho_z, sz, rho_eps, psi;

betta = 0.99;
gama = 2;
alpa = 0.35;
delta = 0.025;
rho_z = 0.95;
sz = 0.01;
rho_eps = 1; % very close (or equal) to 1
psi = 0; % just to provide a initial value

model;

%Law of motion for log TFP
z = rho_z*z(-1) + sz*u;
%epsilon variable
eps = rho_eps*eps(-1);

%budget constraint
0 = exp(ey) - exp(c) - exp(ei);
%output production function
0 = exp(z)*(exp(k(-1))^alpa) - exp(ey);
%law of motion capital
0 = exp(ei) - exp(k) + (1-delta)*exp(k(-1));
%Euler equation:
%Here ((betta/(1 + psi*(1-eps) ) is the overall discount factor. When
%eps=1, it goes back to just betta.
0 = (betta/(1 + psi*(1-eps) ))*exp(-gama*c(+1))*( alpa*exp(ey(+1)-k) + (1-delta)) - exp(-gama*c);

end;

initval;
k = log((alpa/(1/betta - (1-delta)))^(1/(1-alpa)));
z = 0;
eps = 0;
ey = log(exp(k)^alpa);
ei = log(delta*exp(k));  
c = log(exp(ey) - exp(ei)); 
u = 0;
end;

shocks;
var u; stderr 1;
end;

steady;
check;

stoch_simul(order=3,periods=0,irf=0);