The 17 current codes in this folder are:

prepare_model_dsge4.mod: mod file of the two-country DSGE model with bonds and equities.

eval_sss_dsge4_dyn.m: evaluation function to implement the SSS algorithm, with symmetric countries.
eval_sss_dsge4_asym_dyn.m: evaluation function to implement the SSS algorithm, with asymmetric countries.
calib_sss_dsge4_dyn.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4_dyn.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_sss_dsge4_asym_dyn.m: function to approximate the SSS of bonds and equities, with asymmetric countries.
compute_calib_dsge4_dyn.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.

pre_processing_dsge4.m: script to do the pre-processing of the two-country DSGE model.

dsge_ss_y4.m: function to evaluate equation E.15 in Appendix E.2 (to recalculate DSS with asymmetric countries).
euler_errors_dsge4_dyn: function to compute Euler equation errors.
cross_moments_dsge4_dyn: function to compute conditional moments of returns.

solution_dsge4_compstats_integration_dyn.m: script to perform comparative statics (Figure 5).
solution_dsge4_compstats_bonds_dyn.m: script to perform comparative statics (Figure 11 in Appendix E).
solution_dsge4_compstats_asymmetric_dyn.m: script to perform comparative statics (Figure 12 in Appendix E).
solution_dsge4_irfs_dyn.m: script to generate IRFs (Figures 13-15 in Appendix E).
solution_dsge4_simul_dyn.m: script to generate stochastic simulation (Table 3, columns 1-3).
solution_dsge4_simul_asym_dyn.m: script to generate stochastic simulation (Table 3, columns 4 with asymmetric countries).

INSTRUCTIONS: 
- Add the folder "common functions" to the Matlab path.
- Run the script pre_processing_dsge4.m once at the start of each matlab session.