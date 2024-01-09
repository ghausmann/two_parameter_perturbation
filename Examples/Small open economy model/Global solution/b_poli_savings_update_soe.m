function y = b_poli_savings_update_soe(P,b_grid,ui,Pr,b_pol,phi,R)
%This function uses the endogenous grid method, mapping an old saving
%rule to a new one. When both are the same we have the optimal SS rule. 
%Inputs are P (vector of parameters), b_grid and ui (grids for bonds of
%length n_b and income shock values of length n_u), Pr (Probability
%transition matrix of size n_u-by-n_u), b_pol (the old/tomorrow's saving
%rule, a n_u-by-n_b matrix), phi (borrowing limit), and R (gross interest
%rate).
%
% Copyright (C) 2024 Guillermo Hausmann Guil

%Parameter values
beta = P(1);
gamma = P(2);

%lenghts of the asset and income grids
n_b = length(b_grid);
n_u = length(ui);

%Consumption rule tomorrow
c_pol= exp(ui') + b_grid - (1/R )*b_pol;

c_poli = zeros(n_u,n_b);


RHS = R*beta*Pr*((1./c_pol).^gamma);
c_i = (1./RHS).^(1/gamma);
%Compute today's implied initial level of assets.
b_grid1 = c_i + (1/R)*b_grid - exp(ui');

%Interpolate to complete the mapping from the asset grid to the new
%consumption rule
for j=1:n_u
    
    c_poli(j,:) = interp1(b_grid1(j,:), c_i(j,:), b_grid, 'linear', 'extrap');
    
end

%Correct values that violate the constraints:

%Compute the savings rule
b_poli = R*( exp(ui') + b_grid - c_poli);

b_poli(b_poli<=(-phi) )= -phi;
b_poli(b_poli>=(b_grid(n_b))) = b_grid(n_b);

%Return today's policy rule
y = b_poli;