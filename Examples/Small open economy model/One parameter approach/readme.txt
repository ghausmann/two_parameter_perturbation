The 14 current codes in this folder are:

derivatives_soe_o2_p1.m: script to compute up to second order derivatives for the full model
derivatives_soe_o2_p1_only_y.m: script to compute up to second order derivatives for the basic model
derivatives_soe_o3_p1.m: script to compute up to third order derivatives for the full model

sol_soe_o1.m: function to compute the first-order solution of the full model
sol_soe_o1_only_y.m: function to compute the first-order solution of the basic model
sol_soe_o2_only_y.m: function to compute the second-order solution of the basic model
system_soe_o2.m: function to evaluate second-order derivatives of the Euler equation w.r.t. states
system_soe_o3.m: function to evaluate third-order derivatives of the Euler equation w.r.t. states

eval_sss_soe_pac_p1.m: evaluation function to implement the SSS algorithm, full model
eval_sss_soe_pac_p1_only_y.m: evaluation function to implement the SSS algorithm, basic model
calib_sss_soe_pac_p1.m: evaluation function to implement the SSS algorithm in calibration mode, full model
btp_eval.m: function, evaluates the approximated decision rule for bonds

comp_statics_p1_soe.m: script to perform comparative statics
solution_p1_soe.m: script to compute the third-order solution and a simulation


INSTRUCTIONS: 
- Just execute the last two scripts.