function y = c_poli_growth_fast(P,y_grid,k_grid,z,Pr,c_pol)
%This function uses the endogenous gridpoints method to update the
%consumption rule of the planner in the standard growth model. 
%It takes as inputs P (vector of parameters), y_grid (grid for wealth of
%size I), k_grid (grid for capital od size I), z (grid for productivity of
%size S), Pr (probability transition matrix of size SxS), and c_pol
%(current consumption rule of size SxI).
%It returns the updated consumption rule of size SxI.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

betta = P(1);
gama = P(2);
alpa = P(3);
delta = P(4);

I = length(y_grid);
S = length(z);


c_poli = zeros(S,I);
c_pols  = zeros(S,I);

for s=1:S
   
    y_grids = exp(z(s))*k_grid.^(alpa) + (1 -delta)*k_grid;
    c_pols(s,:) = interp1(y_grid,c_pol(s,:),y_grids,'spline','extrap');
end


c_i = ((1/(betta))./(     Pr*( ( alpa*(exp(z')*(k_grid.^(alpa-1))) + 1 -delta     ).*((1./c_pols).^gama))   )     ).^(1/gama);
y_grid1 = k_grid + c_i;


for j=1:S
    
    c_poli(j,:) = interp1(y_grid1(j,:), c_i(j,:), y_grid, 'linear', 'extrap');
   
end

y = c_poli;