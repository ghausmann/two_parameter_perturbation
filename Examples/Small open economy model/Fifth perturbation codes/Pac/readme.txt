The current files in this folder are:

prepare_model_soe_pac.m: script to do the pre-processind of the SOE model (auxiliary model is PAC).

eval_sss_soe_pac.m: function to implement the algorithm described in Section 2.4 of the paper.
calib_sss_soe_pac.m: function to implement the algorithm described in Section 2.4 of the paper, in calibration mode.

compute_sss_soe_pac.m: function to approximate the SSS of bonds.
compute_beta_soe_pac.m: function to calibrate pssi2_star consistnat with a target SSS of bonds.

euler_errors_soe.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_pac.m: script to perform comparative statics.
solution_soe_pac.m: scrpit to compute Euler equation errors and kernel distributions.

INSTRUCTIONS: 
- Add the folders "Perturbation" and "common functions" to the Matlab path.
- Run the script prepare_model_soe_pac.m once. 