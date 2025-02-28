These codes replicate the results from Section 3 (Devereux and Sutherland model) using Levintal's toolkit.

The 13 current codes in this folder are:

prepare_model_ds_5p.m: script, performs the pre-processing of the DS model.

eval_sss_ds_5p.m: evaluation function to implement the SSS algorithms.
eval2_sss_ds_5p.m: evaluation function, returns sum of corrections by epsilon and sigma from Taylor series (see Section 3.6 in paper).

compute_sss_ds_5p.m: function to approximate the SSS of the Home asset.
euler_errors_ds_5p.m: function to compute Euler errors.
errors_routine_5p.m: function, computes expected Euler errors for a given psi.

solution_ds_sss_5p.m: script, replicates Figure 1 in the paper. 
solution_ds_compstats_5p.m: script, replicates Figure 2 in the paper.
solution_ds_irfs_5p.m: script, replicates Figure 3 in the paper.
solution_ds_simul_distribution_5p.m: script, replicates Figure 4 in the paper.
solution_ds_simul_location_5p.m: script, replicates Figure 5 in the paper.
solution_ds_simul_errors_5p.m: script, replicates Table 1 in the paper.

optimal_pssi_ds_5p.m: script, finds psi that maximizes accuracy.

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_ds_5p.m once.

