function y = eval_sss_soe_pac_p1_only_y(b_bar,params)
%SOE with PAC as auxiliary model
%Basic version with only income shocks
%One-parameter perturbation
%This function evaluates the SSS condition (hs2)


%STEP 1: Given b_bar, recalculate the new DSS and implied auxiliary
%parameter values, modifying other function inputs when necessary.
P = params;
%Recalculate DSS
R0 = P(7);
b_bar1 = b_bar;
C01 = 1 + (1-(1/R0))*b_bar1;
P(8) = b_bar;
P(9) = C01;

%STEP 2: Using the new DSS and other function inputs, call an external
%algorithm to compute the matrices H1,..,Hk of derivatives.
%Use closed-form solutions
sol_o1_aux = sol_soe_o1_only_y(P);
sol_o2_aux = sol_soe_o2_only_y(sol_o1_aux,P);

%STEP 3: Using the matrices H1,..,Hk, evaluate the k-order Taylor series of
%h at the DSS, but imposing the model of interest sigma=1.
%Here the shortcut is that if b_bar is the SSS, then hs2=0.
y = sol_o2_aux(1);