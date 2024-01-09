%---------------------------------------------------------------------------
% Small open economy 
% Global solution (endogenous grid method)
%
% Copyright (C) 2024 Guillermo Hausmann Guil
%---------------------------------------------------------------------------

clear;

% Model Parameters
R0 = 1.04; %DSS gross interest rate
betta = 0.9999/R0; %discount factor
gama =  4; %risk aversion
phi_ss =20; %Borrowing limit
rho_y = 0.85; %auto-corr income
ty = 0.02; %conditional std. of income shocks
stdy = ty/((1-rho_y^2)^0.5); %Unconditional std. of income shocks  

P = [betta gama phi_ss];

%Grid-based algorithm parameters
mytol = 1e-14; %Tolerance for convergence
n_b = 1000; %Number of grid points for asset holdings
bmin   = -(phi_ss); % lower bound 
bmax   =50; %upper bound
%Asset grid
b_grid = [bmin bmin + ((1:(n_b-1))/(n_b-1)).^2*(bmax - bmin)];
n_u = 25; %Number of grid points for the idiosyncratic transitory shock
%Use this to approximate the individual iid. shock using the Rouwenhorst's method
[Pr, x_Rouw] = rouwen(rho_y, 0,stdy,n_u);
ui = x_Rouw'; %Individual income shock grid
Pr = Pr'; %Probability transition matrix (rows add up to one)

%Initial guesses
b_pol_guess =  0.9*b_grid + 0.9*ui';  %bond rule guess
b_pol_guess(b_pol_guess<=(-phi_ss)) =  b_grid(1);
b_pol = b_pol_guess;

%Compute the optimal saving rule
dif=1;
while dif > mytol

    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method.
    b_poli = b_poli_savings_update_soe(P,b_grid,ui,Pr,b_pol,phi_ss,R0);
    dif = max(max(abs(b_poli - b_pol)));  %Check covergence
    b_pol = b_poli; %Replace the old rule with the new one

end

%compute the stochastic steady-state
sm = median(1:n_u);
mkpolifun = @(x)interp1(b_grid,b_poli(sm,:),x,'spline','extrap') - x;
bsss = fsolve(mkpolifun,0,optimset('TolFun',1e-15,'Display','off'))





