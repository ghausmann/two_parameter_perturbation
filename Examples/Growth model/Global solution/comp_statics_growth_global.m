% -------------------------------------------------------------------------
% Neoclassical growth model: global solution
% Compute steady-states for different levels of risk, and risk aversion
% -------------------------------------------------------------------------

clear;

%parameters
betta = 0.99;
%gama = 2; %standard calibration
gama = 20; %extreme calibration
alpa = 0.35;
delta = 0.025;
rho_z = 0.95; %persistence
sz = 0.05; %std shocks

%deterministic steady state
Kss = ((1/betta - (1-delta))/alpa)^(-1/(1-alpa));

%Construct grid for capital
I = 400;
Kmin = Kss*0.1;
Kmax = Kss*3;
k_grid = [Kmin Kmin + ((1:(I-1))/(I-1)).^2 * (Kmax - Kmin)];

%Tolerance parameter
mytol = 1e-12;

%Approximate process for exogenous productivity using the Rouwenhorts
%method
S = 15; % size of
sm = median(1:S);
[Pr, zz] = rouwen(rho_z, 0,(sz/((1-rho_z^2)^0.5)),S);
zz = zz'; %grid of productivity values
Pr = Pr'; %Probability transition matrix (rows add up to one)

%Initial guess for law of motion for capital and consumption
k_pol_guess = Kss + 0.9*(k_grid-Kss) + 0.5*zz';
c_pol = exp(zz')*(k_grid.^alpa)  + (1-delta)*k_grid - k_pol_guess;
c_pol(c_pol<=0.01) = 0.01;

%To use the endogenous gridpoints method, define a new state variable
%"wealth"
ymin = exp(zz(1))*(k_grid(1).^alpa)  + (1-delta)*k_grid(1);
ymax = exp(zz(S))*(k_grid(I).^alpa)  + (1-delta)*k_grid(I);
y_grid = [ymin ymin + ((1:(I-1))/(I-1)).^2* (ymax - ymin)];

msz = 0.0025:0.0025:0.03;
lsz = length(msz);

% -------------------------------------------------------------------------
% baseline calibration: gamma = 2
% -------------------------------------------------------------------------

%Vector of parameters
gama1 = 2;
P1 = [betta gama1 alpa delta];
my_ksss1_g = zeros(lsz,1);

for t=1:lsz
    
    [Pr, zz] = rouwen(rho_z, 0,(msz(t)/((1-rho_z^2)^0.5)),S);
    zz = zz'; %grid of productivity values
    Pr = Pr'; %Probability transition matrix (rows add up to one)
    
    %Iterate until convergence
    dif = 1;
    while dif > mytol
        
        c_poli = c_poli_growth_fast(P1,y_grid,k_grid,zz,Pr,c_pol);
        dif = max(max(abs(c_poli - c_pol)));
        c_pol = c_poli;
        
    end
    
    %Put the obtained solution in terms of capital as the state variable
    %(fixing productivity at ss)
    rk_poli = y_grid - c_pol(sm,:);
    myfun = @(y_grid,x)(x.^(alpa)) + (1-delta)*x - y_grid;
    k_val = arrayfun(@(i) fsolve(@(x) myfun(y_grid(i),x),Kss,optimset('Display','off','TolFun',1e-12)),1:numel(y_grid));
    
    %Policy rule for capital
    k_poli = interp1(k_val,rk_poli,k_grid,'spline','extrap'); %global
    
    %Compute stochastic steady-state
    mkpolifun = @(x)interp1(k_grid,k_poli,x,'spline','extrap') - x;
    my_ksss1_g(t) = fsolve(mkpolifun,Kss,optimset('Display','off','TolFun',1e-12));
    
end


% -------------------------------------------------------------------------
% extreme calibration: gamma = 20
% -------------------------------------------------------------------------

gama2 = 20;
P2 = [betta gama2 alpa delta];
my_ksss2_g = zeros(lsz,1);

for t=1:lsz
    
    [Pr, zz] = rouwen(rho_z, 0,(msz(t)/((1-rho_z^2)^0.5)),S);
    zz = zz'; 
    Pr = Pr'; 
    
    dif = 1;
    while dif > mytol
        
        c_poli = c_poli_growth_fast(P2,y_grid,k_grid,zz,Pr,c_pol);
        dif = max(max(abs(c_poli - c_pol)));
        c_pol = c_poli;
        
    end
    
    rk_poli = y_grid - c_pol(sm,:);
    myfun = @(y_grid,x)(x.^(alpa)) + (1-delta)*x - y_grid;
    k_val = arrayfun(@(i) fsolve(@(x) myfun(y_grid(i),x),Kss,optimset('Display','off','TolFun',1e-12)),1:numel(y_grid));
    k_poli = interp1(k_val,rk_poli,k_grid,'spline','extrap'); 
    mkpolifun = @(x)interp1(k_grid,k_poli,x,'spline','extrap') - x;
    my_ksss2_g(t) = fsolve(mkpolifun,Kss,optimset('Display','off','TolFun',1e-12));
    
end

%save results
save('my_comp_sss_global.mat','my_ksss1_g','my_ksss2_g');



