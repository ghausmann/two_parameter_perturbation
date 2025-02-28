These codes replicate the results from Section 6 and Appendix E (multi-asset DSGE model) using Dynare.

The 17 current codes in this folder are:

prepare_model_dsge4.mod: mod file of the two-country DSGE model with bonds and equities.

eval_sss_dsge4.m: evaluation function to implement the SSS algorithm, with symmetric countries.
eval_sss_dsge4_asym.m: evaluation function to implement the SSS algorithm, with asymmetric countries.
calib_sss_dsge4.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_sss_dsge4_asym.m: function to approximate the SSS of bonds and equities, with asymmetric countries.
compute_calib_dsge4.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.

pre_processing_dsge4.m: script to do the pre-processing of the two-country DSGE model.

dsge_ss_y4.m: function to evaluate equation E.26 in Appendix E.2 (to recalculate DSS with asymmetric countries).
euler_errors_dsge4: function to compute Euler equation errors.
cross_moments_dsge4: function to compute conditional moments of returns.

solution_dsge4_compstats_integration.m: script to perform comparative statics (Figure 6).
solution_dsge4_compstats_bonds.m: script to perform comparative statics (Figure 6 in Appendix E).
solution_dsge4_compstats_asymmetric.m: script to perform comparative statics (Figure 7 in Appendix E).
solution_dsge4_irfs.m: script to generate IRFs (Figures 8-10 in Appendix E).
solution_dsge4_simul.m: script to generate stochastic simulation (Table 3, columns 1-3).
solution_dsge4_simul_asym.m: script to generate stochastic simulation (Table 3, column 4 with asymmetric countries).

INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_dsge4.m once at the start of each Matlab session.