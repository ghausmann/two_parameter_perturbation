function y = calib_sss_soe_pac_p1(psi2,b_bar,params)
%SOE with PAC as auxiliary model
%Full version with income and interest rate shocks
%One-parameter perturbation
%This function evaluates the SSS condition (hs2) in calibration mode

%STEP 1: Given psi2 and b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.
P = params;
%Assign new value to psi2
P(11) = psi2;
%Recalculate DSS
R0 = P(7);
b_bar1 = b_bar;
C01 = 1 + (1-(1/R0))*b_bar1;
P(8) = b_bar;
P(9) = C01;

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices H1,..,Hk of derivatives.
%Use closed-form solution for first-order, and non-linear solver for second
sol_o1_aux = sol_soe_o1(P);
my_soe_system = @(x)system_soe_o2(x,sol_o1_aux,P);
sol_o2 = fsolve(my_soe_system,zeros(1,7),optimset('TolFun',1e-12,'Display','off'));

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest sigma=1.
%Here the shortcut is that if b_bar is the SSS, then hs2=0.
y = sol_o2(1);