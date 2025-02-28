These codes replicate results from Appendix D.3 (SOE model) using Dynare.
The auxiliary model uses Uzawa preferences.

The 11 current codes in this folder are:

soe_uzawa_basic.mod: mod file of the SOE (Uzawa as auxiliary model), with the calibration for comparative statics.
soe_uzawa_full.mod: mod file of the SOE (Uzawa as auxiliary model), with the calibration using mexican data

eval_sss_soe_uzawa_dynare.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_uzawa_dynare.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_soe_uzawa_dynare.m: function to approximate the SSS of bonds.
compute_beta_soe_uzawa_dynare.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

pre_processing_soe_uzawa_basic.m: script to do the pre-processing of the SOE model (calibration for comparative statics).
pre_processing_soe_uzawa_full.m: script to do the pre-processing of the SOE model (calibration using mexican data).

euler_errors_soe_dynare.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_uzawa_dynare.m: script to perform comparative statics (Figure 3 from Appendix D).
solution_soe_uzawa_dynare.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D).


INSTRUCTIONS: 
- Add the folder "common functions" to the Matlab path.
- Run the scripts pre_processing_soe_uzawa_basic.m and pre_processing_soe_uzawa_full.m once at the start of each matlab session.