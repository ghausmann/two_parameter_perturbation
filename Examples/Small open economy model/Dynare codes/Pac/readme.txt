The current files in this folder are:

soe_pac_basic.mod: mod file of the SOE (PAC as auxiliary model), with the calibration for comparative statics.
soe_pac_full.mod: mod file of the SOE (PAC as auxiliary model), with the calibration using mexican data

eval_sss_soe_pac_dynare.m: function to implement the algorithm described in Section 2.4 of the paper.
calib_sss_soe_pac_dynare.m: function to implement the algorithm described in Section 2.4 of the paper, in calibration mode.

compute_sss_soe_pac_dynare.m: function to approximate the SSS of bonds.
compute_beta_soe_pac_dynare.m: function to calibrate pssi2_star consistnat with a target SSS of bonds.

pre_processing_soe_pac_basic.m: script to do the pre-processing of the SOE model (calibration for comparative statics).
pre_processing_soe_pac_full.m: script to do the pre-processing of the SOE model (calibration using mexican data).

euler_errors_soe_dynare.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_pac_dynare.m: script to perform comparative statics.
solution_soe_pac_dynare.m: scrpit to compute Euler equation errors and kernel distributions.


INSTRUCTIONS: 
- Add the folder "common functions" to the Matlab path.
- Run the scripts pre_processing_soe_pac_basic.m and pre_processing_soe_pac_full.m once at the start of each matlab session.