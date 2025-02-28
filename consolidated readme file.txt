Consolidated (folder by folder) readme file.

It includes approximate computing times for scripts using my personal laptop (four years old, with a 2.10GHz AMD Ryzen 5 5500U processor with six cores).

%-------------------------------------------------
1) Devereux-Sutherland model\ds model dynare
%-------------------------------------------------

These codes replicate the results from Section 3 (Devereux and Sutherland model) using Dynare.

The 12 current codes in this folder are:

prepare_model_ds.mod: Dynare MOD file
pre_processing_ds.m: script, performs the pre-processing of the DS model. (3 seconds)

eval_sss_ds.m: evaluation function to implement the SSS algorithms.
eval2_sss_ds.m: evaluation function, returns sum of corrections by epsilon and sigma from Taylor series (see Section 3.6 in paper).

compute_sss_ds.m: function to approximate the SSS of the Home asset.
euler_errors_ds.m: function to compute Euler errors.

solution_ds_sss.m: script, replicates Figure 1 in the paper. (0.41 seconds)
solution_ds_compstats.m: script, replicates Figure 2 in the paper. (2 seconds)
solution_ds_irfs.m: script, replicates Figure 3 in the paper. (2 seconds)
solution_ds_simul_distribution.m: script, replicates Figure 4 in the paper. (36 seconds)
solution_ds_simul_location.m: script, replicates Figure 5 in the paper. (50 seconds)
solution_ds_simul_errors.m: script, replicates Table 1 in the paper (up to third-order). (3 seconds)


INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_ds.m once at the start of each Matlab session.

%-------------------------------------------------
2) Devereux-Sutherland model\ds model 5p
%-------------------------------------------------

These codes replicate the results from Section 3 (Devereux and Sutherland model) using Levintal's toolkit.

The 13 current codes in this folder are:

prepare_model_ds_5p.m: script, performs the pre-processing of the DS model. (3 seconds)

eval_sss_ds_5p.m: evaluation function to implement the SSS algorithms.
eval2_sss_ds_5p.m: evaluation function, returns sum of corrections by epsilon and sigma from Taylor series (see Section 3.6 in paper).

compute_sss_ds_5p.m: function to approximate the SSS of the Home asset.
euler_errors_ds_5p.m: function to compute Euler errors.
errors_routine_5p.m: function, computes expected Euler errors for a given psi.

solution_ds_sss_5p.m: script, replicates Figure 1 in the paper. (0.92 seconds)
solution_ds_compstats_5p.m: script, replicates Figure 2 in the paper. (2 seconds)
solution_ds_irfs_5p.m: script, replicates Figure 3 in the paper. (0.97 seconds)
solution_ds_simul_distribution_5p.m: script, replicates Figure 4 in the paper. (56 seconds)
solution_ds_simul_location_5p.m: script, replicates Figure 5 in the paper. (7 seconds)
solution_ds_simul_errors_5p.m: script, replicates Table 1 in the paper. (3 seconds the third-order solution, 26 seconds the fourth-order, 200 seconds the fifth-order)

optimal_pssi_ds_5p.m: script, finds psi that maximizes accuracy. (4 seconds)

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_ds_5p.m once.


%-------------------------------------------------
3) Devereux-Sutherland model\ds model wealth 5p
%-------------------------------------------------

These codes replicate the results from Appendix A (Devereux and Sutherland model) using Levintal's toolkit.

The 5 current codes in this folder are:

prepare_model_ds_wealth_5p.m: script, performs the pre-processing of the DS model, using the original notation by DS. (2 seconds)

eval_sss_ds_wealth_5p.m: evaluation function to implement the SSS algorithm.
compute_sss_wealth_ds_5p.m: function to approximate the SSS

solution_ds_compstats_wealth_5p.m: script, replicates Figure 1, panel (a) in Appendix A. (2 seconds)
solution_ds_irfs_wealth_5p.m: script, replicates Figure 1, panel (b) in Appendix A. (0.80 seconds)

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_ds_wealth_5p.m once.


%--------------------------------------------------------
4) multi-asset DSGE model\cg model dynare\Equities only
%--------------------------------------------------------

These codes replicate results from Appendix E (multi-asset DSGE model) using Dynare.

The 5 current codes in this folder are:

prepare_model_dsge2.mod: mod file of the two-country DSGE model with equities only.
eval_sss_dsge2.m: evaluation function to implement the SSS algorithm.
compute_sss_dsge2.m: function, approximates the SSS of equities.

pre_processing_dsge2.m: script, pre-processing of the two-country DSGE model. (2 seconds)

solution_dsge2_compstatsn.m: script, perform comparative statics (Figure 6 in Appendix E). (1 second)

INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_dsge2.m once at the start of each Matlab session.


%--------------------------------------------------------------
5) multi-asset DSGE model\cg model dynare\Bonds and equities
%--------------------------------------------------------------

These codes replicate the results from Section 6 and Appendix E (multi-asset DSGE model) using Dynare.

The 17 current codes in this folder are:

prepare_model_dsge4.mod: mod file of the two-country DSGE model with bonds and equities. 

eval_sss_dsge4.m: evaluation function to implement the SSS algorithm, with symmetric countries.
eval_sss_dsge4_asym.m: evaluation function to implement the SSS algorithm, with asymmetric countries.
calib_sss_dsge4.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_sss_dsge4_asym.m: function to approximate the SSS of bonds and equities, with asymmetric countries.
compute_calib_dsge4.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.

pre_processing_dsge4.m: script to do the pre-processing of the two-country DSGE model. (2 seconds)

dsge_ss_y4.m: function to evaluate equation E.26 in Appendix E.2 (to recalculate DSS with asymmetric countries).
euler_errors_dsge4: function to compute Euler equation errors.
cross_moments_dsge4: function to compute conditional moments of returns.

solution_dsge4_compstats_integration.m: script to perform comparative statics (Figure 6). (3 seconds)
solution_dsge4_compstats_bonds.m: script to perform comparative statics (Figure 6 in Appendix E). (6 seconds)
solution_dsge4_compstats_asymmetric.m: script to perform comparative statics (Figure 7 in Appendix E). (15 seconds)
solution_dsge4_irfs.m: script to generate IRFs (Figures 8-10 in Appendix E). (2 seconds each figure)
solution_dsge4_simul.m: script to generate stochastic simulation (Table 3, columns 1-3). (51 seconds each column)
solution_dsge4_simul_asym.m: script to generate stochastic simulation (Table 3, column 4 with asymmetric countries). (55 seconds)

INSTRUCTIONS: 
- Add the folder "common functions" and Dynare's (5x version) subfolder "matlab" (folder and subfolders) to the Matlab path.
- Run the script pre_processing_dsge4.m once at the start of each Matlab session.


%--------------------------------------------------------------
6) multi-asset DSGE model\cg model 5p\Bonds and equities
%--------------------------------------------------------------

These codes replicate the results from Section 6 and Appendix E (multi-asset DSGE model) using Levintal's toolkit

The 11 current codes in this folder are:

prepare_model_dsge4_5p.m: script, performs the pre-processing of the two-country DSGE model with bonds and equities. (4 seconds)
prepare_model_dsge4_eps1_5p.m: same as prepare_model_dsge4_5p.m, but the perturbation parameter epsilon is NOT squared. (4 seconds)  

eval_sss_dsge4_5p.m: evaluation function to implement the SSS algorithm, with symmetric countries. 
calib_sss_dsge4_5p.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_dsge4_5p.m: function to approximate the SSS of bonds and equities, with symmetric countries.
compute_calib_dsge4_5p.m: function to approximate the SSS of bonds, and std. of preference shocks consistent with a target value of equities.
euler_errors_dsge4_5p.m: function to compute Euler errors for higher-order solutions
euler_errors_dsge4_linear_5p.m: function to compute Euler errors for linear solutions (both first-order and risk-linear)

solution_dsge4_compstats_integration_5p.m: script to perform comparative statics (Figure 6 in paper). (3 seconds)
solution_dsge4_simul_5p.m: script to generate the stochastic simulation (Table 3, columns 1-3). (59 seconds each column)
solution_dsge4_robustness_5p.m: script to generate the robustness checks (Table 4 in Appendix E). (59 seconds columns 1,5,6; 4 seconds column 2,4; 15 seconds column 3)

INSTRUCTIONS: 
- Add the folders "Perturbation" (folder and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_dsge4_5p.m once.


%--------------------------------------------------------------
7) Examples\Growth model\Dynare codes
%--------------------------------------------------------------

These codes replicate results from Appendix D.2 (Growth model) using Dynare.

The 5 current codes in this folder are:

prepare_model_growth.mod: mod file describing the model.
pre_processing_growth_dyn.m: script, pre-processing of the growth model. (2 seconds)

eval_sss_growth_dyn.m: evaluation function, implements the SSS algorithm.
compute_sss_growth_dyn.m: function, approximates the SSS of capital using two-parameter perturbation.

simulation_growth_dyn.m: script, simulates the model (Table 1 in Appendix D.2). (20 seconds)

And the current data files (from the global solution) are:
my_k_poli.mat

INSTRUCTIONS: 
- Add the main folder "common functions" to the Matlab path.
  If you plan to change the default scenario and compare results with the global solution, add the subfolder "Growth model\Global solution" to the Matlab path too.
- Run the script pre_processing_growth_dyn.m once at the start of each MATLAB session.


%--------------------------------------------------------------
8) Examples\Growth model\Fifth perturbation codes
%--------------------------------------------------------------

These codes replicate results from Appendix D.2 (Growth model) using Levintal's toolkit.

The 6 current codes in this folder are:

prepare_model_growth.m: script, pre-processing of the growth model. (0.91 seconds)

eval_sss_growth.m: evaluation function, implements the SSS algorithm.
compute_sss_growth.m: function, approximates the SSS of capital using two-parameter perturbation.
iter_sss_growth: function, approximates SSS of capital using standard perturbation.

comp_statics_growth.m: script, performs comparative statics (Figure 2 in Appendix D.2). (1 second)
simulation_growth.m: script, simulates the model (Table 1 in Appendix D.2). (1 second)

And the current data files (from the global solution) are:
my_comp_sss_global.mat
my_k_poli.mat

INSTRUCTIONS: 
- Add the main folders "Perturbation" (folders and subfolders) and "common functions" to the Matlab path. 
  If you plan to change the default scenario and compare results with the global solution, add the subfolder "Growth model\Global solution" to the Matlab path too. 
- Run the script prepare_model_growth.m once. 


%--------------------------------------------------------------
9) Examples\Growth model\Global solution
%--------------------------------------------------------------

These codes replicate results from Appendix D.2 (Growth model) using a global solution method.

The current codes in this folder are:

c_poli_growth_fast.m: function, implements the endogenous grid method.
solution_growth_global.m: script, solves the growth model. (8 seconds)
comp_statics_growth_global.m: script, performs comparative statics. (53 seconds)

INSTRUCTIONS: 
- Add the main folder "common functions" to the Matlab path.


%--------------------------------------------------------------
10) Examples\Small open economy model\Dynare codes\Pac
%--------------------------------------------------------------

These codes replicate results from Appendix D.3 (SOE model) using Dynare.
The auxiliary model uses PAC.

The 12 current codes in this folder are:

soe_pac_basic.mod: mod file of the SOE (PAC as auxiliary model), with the calibration for comparative statics.
soe_pac_full.mod: mod file of the SOE (PAC as auxiliary model), with the calibration using mexican data

eval_sss_soe_pac_dynare.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_pac_dynare.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_soe_pac_dynare.m: function to approximate the SSS of bonds.
compute_beta_soe_pac_dynare.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

pre_processing_soe_pac_basic.m: script to do the pre-processing of the SOE model (calibration for comparative statics). (2 seconds)
pre_processing_soe_pac_full.m: script to do the pre-processing of the SOE model (calibration using mexican data). (2 seconds)

euler_errors_soe_dynare.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_pac_dynare.m: script to perform comparative statics. (Figure 3 from Appendix D). (3 seconds)
solution_soe_pac_dynare.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D). (8 seconds)
solution_soe_pac_compare.m: script to compare two-parameter perturbation with a single perturbation parameter approach. (2 seconds)


INSTRUCTIONS: 
- Add the folder "common functions" to the Matlab path.
- Run the scripts pre_processing_soe_pac_basic.m and pre_processing_soe_pac_full.m once at the start of each Matlab session.


%--------------------------------------------------------------
11) Examples\Small open economy model\Dynare codes\Uzawa
%--------------------------------------------------------------

These codes replicate results from Appendix D.3 (SOE model) using Dynare.
The auxiliary model uses Uzawa preferences.

The 11 current codes in this folder are:

soe_uzawa_basic.mod: mod file of the SOE (Uzawa as auxiliary model), with the calibration for comparative statics.
soe_uzawa_full.mod: mod file of the SOE (Uzawa as auxiliary model), with the calibration using mexican data

eval_sss_soe_uzawa_dynare.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_uzawa_dynare.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_soe_uzawa_dynare.m: function to approximate the SSS of bonds.
compute_beta_soe_uzawa_dynare.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

pre_processing_soe_uzawa_basic.m: script to do the pre-processing of the SOE model (calibration for comparative statics). (2 seconds)
pre_processing_soe_uzawa_full.m: script to do the pre-processing of the SOE model (calibration using mexican data). (2 seconds)

euler_errors_soe_dynare.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_uzawa_dynare.m: script to perform comparative statics (Figure 3 from Appendix D). (3 seconds)
solution_soe_uzawa_dynare.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D). (7 seconds)


INSTRUCTIONS: 
- Add the folder "common functions" to the Matlab path.
- Run the scripts pre_processing_soe_uzawa_basic.m and pre_processing_soe_uzawa_full.m once at the start of each matlab session.


%------------------------------------------------------------------
12) Examples\Small open economy model\Fifth perturbation codes\Pac
%------------------------------------------------------------------

These codes replicate results from Appendix D.3 (SOE model) using Levintal's toolkit.
The auxiliary model uses PAC.

The 8 current codes in this folder are:

prepare_model_soe_pac.m: script, performs the pre-processing of the SOE model (auxiliary model is PAC). (1 second)

eval_sss_soe_pac.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_pac.m: evaluation function to implement the SSS algorithm, in calibration mode.

compute_sss_soe_pac.m: function to approximate the SSS of bonds.
compute_beta_soe_pac.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

euler_errors_soe.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_pac.m: script to perform comparative statics (Figure 3 from Appendix D). (3 seconds)
solution_soe_pac.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D). (11 seconds)

INSTRUCTIONS: 
- Add the folders "Perturbation" (folders and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_soe_pac.m once. 

%---------------------------------------------------------------------
13) Examples\Small open economy model\Fifth perturbation codes\Uzawa
%---------------------------------------------------------------------

These codes replicate results from Appendix D.3 (SOE model) using Levintal's toolkit.
The auxiliary model uses Uzawa preferences.

The 8 current codes in this folder are:

prepare_model_soe_uzawa.m: script, performs the pre-processing of the SOE model (auxiliary model is Uzawa). (1 second)

eval_sss_soe_uzawa.m: evaluation function to implement the SSS algorithm.
calib_sss_soe_uzawa.m: evaluation function to implement the SSS algorithm in calibration mode.

compute_sss_soe_uzawa.m: function to approximate the SSS of bonds.
compute_beta_soe_uzawa.m: function to calibrate pssi2_star consistent with a target SSS of bonds.

euler_errors_soe.m: function to compute Euler equation errors of the SOE.

comp_statics_soe_uzawa.m: script to perform comparative statics (Figure 3 from Appendix D). (4 seconds)
solution_soe_uzawa.m: script to compute Euler equation errors and kernel distributions (Figure 4 and Table 3 from Appendix D). (11 seconds)

INSTRUCTIONS: 
- Add the folders "Perturbation" (folders and subfolders) and "common functions" to the Matlab path.
- Run the script prepare_model_soe_uzawa.m once. 


%---------------------------------------------------------------------
14) Examples\Small open economy model\Global solution
%---------------------------------------------------------------------

These codes replicate results from Appendix D.3 (SOE model) using a global solution method.

The 3 current codes in this folder are:

b_poli_savings_update_soe.m: function to implement endogenous grid method.
savings_solution_soe.m: script to solve the SOE model with the endogenous grid method. (3 seconds)
comp_statics_soe_global.m: script to perform comparative statics of the SOE model with the endogenous grid method. (6 seconds)


INSTRUCTIONS: 
- Add the folders "common functions" to the Matlab path.


%---------------------------------------------------------------------
15) Examples\Small open economy model\One parameter approach
%---------------------------------------------------------------------

These codes implement a one-parameter approach to solve the SOE of Appendix D.3,
just to show that it delivers identical results to those obtained with two-parameter perturbation.

The 14 current codes in this folder are:

derivatives_soe_o2_p1.m: script to compute up to second order derivatives for the full model
derivatives_soe_o2_p1_only_y.m: script to compute up to second order derivatives for the basic model
derivatives_soe_o3_p1.m: script to compute up to third order derivatives for the full model

sol_soe_o1.m: function to compute the first-order solution of the full model
sol_soe_o1_only_y.m: function to compute the first-order solution of the basic model
sol_soe_o2_only_y.m: function to compute the second-order solution of the basic model
system_soe_o2.m: function to evaluate second-order derivatives of the Euler equation w.r.t. states
system_soe_o3.m: function to evaluate third-order derivatives of the Euler equation w.r.t. states

eval_sss_soe_pac_p1.m: evaluation function to implement the SSS algorithm, full model
eval_sss_soe_pac_p1_only_y.m: evaluation function to implement the SSS algorithm, basic model
calib_sss_soe_pac_p1.m: evaluation function to implement the SSS algorithm in calibration mode, full model
btp_eval.m: function, evaluates the approximated decision rule for bonds

comp_statics_p1_soe.m: script to perform comparative statics (0.36 seconds)
solution_p1_soe.m: script to compute the third-order solution and a simulation (0.08 seconds)


INSTRUCTIONS: 
- Just execute the last two scripts.


%---------------------------------------------------------------------
16) Examples\Small open economy model\One bond, two-country model
%---------------------------------------------------------------------

Codes to solve the one-bond, two-country model of Appendix D.4, using Levintal's toolkit.

The 7 current codes in this folder are:

prepare_mode_bond2.m: script, pre-processing of the model (1 second)
solution_bond2.m: script, replicates Figure 5 in Appendix D.4 (30 seconds)

eval_sss_bond2.m: function, evaluates SSS condition
calib_sss_bond2.m: function, evaluates SSS condition in calibration mode
compute_sss_bond2.m: function, implements SSS algorithm
compute_calib_bond2.m: function, implements SSS algorithm in calibration mode
euler_errors_bond2.m: function, computes Euler errors

INSTRUCTIONS: add the folders "Perturbation" and "common functions" to MATLAB's search path, and run "prepare_mode_bond2.m" once.





