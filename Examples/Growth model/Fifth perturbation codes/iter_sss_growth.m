function y = iter_sss_growth(sz,model,params,M,approx)
% This function iterates over the states without shocks to find the
% stochastic steady-state of the growth model solved by standard
% perturbation.

% Make eta
eta=[0;sz;0];
params(6) = sz;

betta = params(1);
alpa = params(3);
delta = params(4);

% deterministic steady-state
Kss_0 = (alpa/(1/betta - (1-delta)))^(1/(1-alpa));
k1 = Kss_0;
kss = log(k1);
eyss = log(exp(kss)^alpa);
eiss = log(delta*exp(kss));
css = log(exp(eyss) - exp(eiss));
zss = 0;
epsss = 0;

% collect steady-state
nxss=[kss;zss;epsss];
nyss=[css;eiss;eyss];

%With standard perturbation, no role for any auxiliary parameter.
psi = 0;
params(8) = psi;

% choose algorithm
algo='gensylv';
% compute the perturbation solution
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

x0=nxss; % start at the steady state
xt = x0;

% iterate without shocks until convergence to the stochastic steady-state
dif=1;
while dif>1e-15

    xt1=dr_ht(derivs1,nxss,approx,(xt-nxss));
    dif=(xt1(1)-xt(1))^2;
    xt=xt1;
end

% return the implied stochatic steady-state
y = exp(xt1(1));