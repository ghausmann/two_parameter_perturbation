function y = errors_routine_5p(x_pssi,model,params,eta,M,eps_ind,approx,epsi_nodes,weight_nodes,P)
%This function takes as an input a candidate pssi (the parameter scaling
%PAC-like modifications in the Devereux-Sutherland model), and uses it to
%compute expected Euler errors at the SSS (see section 4.2 in the paper).
%It can be used by an optimization routine to maximize accuracy.

z0 = P(5);
params(4) = exp(x_pssi); %in exp to make sure we look for a positive value of pssi

%First, compute SSS
a_sss = compute_sss_ds_5p(model,params,M,eps_ind,0);

%Recalculate DSS of auxiliary model
ah1 = a_sss;
af1 = -a_sss;
params(1) = a_sss;
%params(4) = pssi_v(n);

%Deterministic steady-state
nxss=[ah1;af1;zeros(6,1)]; %for states
nyss = [1;1;z0;z0;1;1]; %for controls
algo='vectorize';
%Compute derivatives of policy rules
derivs1=solve_dsge(model,params,M,eta,nxss,nyss,approx,algo);

%--------------------------------------------------------------------------
% Simulation of the model of interest
%--------------------------------------------------------------------------
x0=nxss; % start at the steady state
x0(eps_ind) = 1;
derivs1.hx(eps_ind,eps_ind) = 1; % evaluate at the model of interest (epsilon=1)

%Metric used is expected portfolio euler errors conditional on being at the
%SSS
[NN ,~] = size(epsi_nodes);
my_errors =zeros(2,NN);
my_x0 = x0 + eta*epsi_nodes';


for nn=1:NN
    
    my_errors(:,nn) = (euler_errors_ds_5p(P,nxss,nyss,my_x0(:,nn),epsi_nodes,weight_nodes,derivs1,eta,approx));
    
end

my_weighted_errors = weight_nodes'*my_errors';

%This output makes it easier for the optimization routine to find the
%minimum.
y =  1/abs(log10(my_weighted_errors(2)));


