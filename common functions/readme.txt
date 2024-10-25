The 15 current codes in this folder are:

dr_gt.m: function, evaluates decision rules for control variables of a DSGE model solved with the algorithm by Levintal (2017).
dr_ht.m: function, evaluates decision rules for state variables of a DSGE model solved with the algorithm by Levintal (2017).
dr_yt.m: function, evaluates decision rules for all variables of a DSGE model solved with Dynare.
mykron1.m: function, computes kronecker products.

my_jac_gx.m: function, computes the jacobian of k-order decision rules for control variables of a DSGE model solved with the algorithm by Levintal (2017).
my_jac_hx.m: function, computes the jacobian of k-order decision rules for state variables of a DSGE model solved with the algorithm by Levintal (2017).
my_jac_yt_dyn.m: function, computes the jacobian of 3-order decision rules for the endogenous variables of a DSGE model solved with Dynare.

simul_linear_risk.m: function, implements a simulation using linear risk-corrected rules obtained with the algorithm by Levintal (2017).
simul_linear_risk_dyn.m: function, implements a simulation using linear risk-corrected rules obtained with Dynare.
simul_mod_pruning3.m: function, implements a pruned simulation treating perturbation objects as parameters.
simult_mod: slightly modified version of the original Dynare function "simult_.m". Implements a pruned simulation treating the perturbation object sigma as a parameter.

Monomials_1.m: external function by Judd, Maliar and Maliar (2011), constructs monomials.
Monomials_2.m: external function by Judd, Maliar and Maliar (2011), constructs monomials.
rouwen.m: external function, implements the Rouwenhorst method.
hitm_s.m: external function, simulates a time-series using a Markov chain.
