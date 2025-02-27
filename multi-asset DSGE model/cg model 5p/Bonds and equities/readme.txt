The 11 current codes in this folder are:

prepare_model_dsge4_5p.m: script, performs the pre-processing of the two-country DSGE model with bonds and equities.
prepare_model_dsge4_eps1_5p.m: same as prepare_model_dsge4_5p.m, but the perturbation parameter epsilon is NOT squared.

eval_sss_dsge4_5p.m: evaluation function to implement the SSS algorithm, with symmetric countries.
calib_sss_dsge4_5p.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4_5p.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_calib_dsge4_5p.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.
euler_errors_dsge4_5p.m: function to compute Euler errors for higher-order solutions
euler_errors_dsge4_linear_5p.m: function to compute Euler errors for linear solutions (both first-order and risk-linear)

solution_dsge4_compstats_integration_5p.m: script to perform comparative statics (Figure 6 in paper).
solution_dsge4_simul_5p.m: script to generate the stochastic simulation (Table 3, columns 1-3).
solution_dsge4_robustness_5p.m: script to generate the robustness checks (Table 4 in Appendix E).

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_dsge4_5p.m once.