These codes replicate results from Appendix E (multi-asset DSGE model) using Dynare.

The 5 current codes in this folder are:

prepare_model_dsge2.mod: mod file of the two-country DSGE model with equities only.
eval_sss_dsge2.m: evaluation function to implement the SSS algorithm.
compute_sss_dsge2.m: function, approximates the SSS of equities.

pre_processing_dsge2.m: script, pre-processing of the two-country DSGE model.

solution_dsge2_compstatsn.m: script, perform comparative statics (Figure 6 in Appendix E).

INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_dsge2.m once at the start of each Matlab session.