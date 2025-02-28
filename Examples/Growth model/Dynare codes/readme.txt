These codes replicate results from Appendix D.2 (Growth model) using Dynare.

The 5 current codes in this folder are:

prepare_model_growth.mod: mod file describing the model.
pre_processing_growth_dyn.m: script, pre-processing of the growth model.

eval_sss_growth_dyn.m: evaluation function, implements the SSS algorithm.
compute_sss_growth_dyn.m: function, approximates the SSS of capital using two-parameter perturbation.

simulation_growth_dyn.m: script, simulates the model (Table 1 in Appendix D.2).

And the current data files (from the global solution) are:
my_k_poli.mat

INSTRUCTIONS: 
- Add the main folder "common functions" to the Matlab path.
  If you plan to change the default scenario and compare results with the global solution, add the subfolder "Growth model\Global solution" to the Matlab path too.
- Run the script pre_processing_growth_dyn.m once at the start of each MATLAB session.
