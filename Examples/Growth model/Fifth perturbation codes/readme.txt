These codes replicate results from Appendix D.2 (Growth model) using Levintal's toolkit.

The 6 current codes in this folder are:

prepare_model_growth.m: script, pre-processing of the growth model.

eval_sss_growth.m: evaluation function, implements the SSS algorithm.
compute_sss_growth.m: function, approximates the SSS of capital using two-parameter perturbation.
iter_sss_growth: function, approximates SSS of capital using standard perturbation.

comp_statics_growth.m: script, performs comparative statics (Figure 2 in Appendix D.2).
simulation_growth.m: script, simulates the model (Table 1 in Appendix D.2).

And the current data files (from the global solution) are:
my_comp_sss_global.mat
my_k_poli.mat

INSTRUCTIONS: 
- Add the main folders "Perturbation" (folders and subfolders) and "common functions" to the Matlab path. 
  If you plan to change the default scenario and compare results with the global solution, add the subfolder "Growth model\Global solution" to the Matlab path too. 
- Run the script prepare_model_growth.m once. 