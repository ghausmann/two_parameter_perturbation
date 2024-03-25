% -------------------------------------------------------------------------
% Neoclassical growth model: global solution
% -------------------------------------------------------------------------

clear;

% Choose the calibration. Zero for standard, and one for extreme.
% Default is extreme calibration.
calibration = 1;

%parameters
betta = 0.99;
alpa = 0.35;
delta = 0.025;
rho_z = 0.95;

if calibration==0
    gama = 2;
    sz = 0.01;
elseif calibration==1
    gama = 20;
    sz = 0.025;
end

%Vector of parameters
P = [betta gama alpa delta];

%deterministic steady state
Kss = ((1/betta - (1-delta))/alpa)^(-1/(1-alpa))

%Construct grid for capital 
I = 400;
Kmin = Kss*0.1;
Kmax = Kss*3;
k_grid = [Kmin Kmin + ((1:(I-1))/(I-1)).^2 * (Kmax - Kmin)];

%Tolerance parameter
mytol = 1e-12;

%Approximate process for exogenous productivity using the Rouwenhort
%method
S = 15; 
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

%Iterate until convergence
dif = 1;
while dif > mytol

    c_poli = c_poli_growth_fast(P,y_grid,k_grid,zz,Pr,c_pol);
    dif = max(max(abs(c_poli - c_pol)));
    c_pol = c_poli;

end

k_poli = zeros(S,I);

%Put the obtained solution in terms of capital as the state variable
for s=1:S
   
    rk_poli = y_grid - c_pol(s,:);
    myfun = @(y_grid,x)exp(zz(s))*(x.^(alpa)) + (1-delta)*x - y_grid;
    k_val = arrayfun(@(i) fsolve(@(x) myfun(y_grid(i),x),Kss,optimset('Display','off')),1:numel(y_grid));
    k_poli(s,:) = interp1(k_val,rk_poli,k_grid,'spline','extrap'); 

end

%figure;
%plot(k_grid,k_grid,'g','LineWidth',1);hold on ;plot(k_grid,k_poli(sm,:),'b','LineWidth',1);
% title('Policy rule when z=0');

% compute stochastic steady-state
mkpolifun = @(x)interp1(k_grid,k_poli(sm,:),x,'spline','extrap') - x;
my_ss = fsolve(mkpolifun,Kss,optimset('TolFun',1e-15))

%save results
save('my_k_poli.mat','k_poli','k_grid','zz','Pr');


