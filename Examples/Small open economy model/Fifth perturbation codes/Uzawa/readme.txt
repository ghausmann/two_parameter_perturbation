The 8 current codes in this folder are:

prepare_model_soe_uzawa.m: script, performs the pre-processing of the SOE model (auxiliary model is Uzawa).

eval_sss_soe_uzawa.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_uzawa.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_soe_uzawa.m: function to approximate the SSS of bonds.
compute_beta_soe_uzawa.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

euler_errors_soe.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_uzawa.m: script to perform comparative statics (Figure 3 from Appendix D).
solution_soe_uzawa.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D).

INSTRUCTIONS: 
- Add the folders "Perturbation" (folders and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_soe_uzawa.m once. 