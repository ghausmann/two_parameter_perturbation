The 8 current codes in this folder are:

prepare_model_dsge4.m: script, performs the pre-processing of the two-country DSGE model with bonds and equities.

eval_sss_dsge4.m: evaluation function to implement the SSS algorithm, with symmetric countries.
calib_sss_dsge4.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_calib_dsge4.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.
euler_errors_dsge4.m: function to compute Euler errors

solution_dsge4_compstats_integration.m: script to perform comparative statics (Figure 5 in paper).
solution_dsge4_simul.m: script to generate the stochastic simulation (Table 3, columns 1-3).

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_dsge4.m once.