These codes replicate the results from Section 3 (Devereux and Sutherland model) using Dynare.

The 12 current codes in this folder are:

prepare_model_ds.mod: Dynare MOD file
pre_processing_ds.m: script, performs the pre-processing of the DS model.

eval_sss_ds.m: evaluation function to implement the SSS algorithms.
eval2_sss_ds.m: evaluation function, returns sum of corrections by epsilon and sigma from Taylor series (see Section 3.6 in paper).

compute_sss_ds.m: function to approximate the SSS of the Home asset.
euler_errors_ds.m: function to compute Euler errors.

solution_ds_sss.m: script, replicates Figure 1 in the paper. 
solution_ds_compstats.m: script, replicates Figure 2 in the paper.
solution_ds_irfs.m: script, replicates Figure 3 in the paper.
solution_ds_simul_distribution.m: script, replicates Figure 4 in the paper.
solution_ds_simul_location.m: script, replicates Figure 5 in the paper.
solution_ds_simul_errors.m: script, replicates Table 1 in the paper (up to third-order).


INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_ds.m once at the start of each Matlab session.