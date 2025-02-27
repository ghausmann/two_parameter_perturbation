% Simple script to do the pre-processing of the SOE model (calibration with
% Mexican data).
% With Dynare, you need to execute this code once at the start of each
% session.
% Also, add the folder "common functions" to the MATLAB path.

clear;
addpath('C:\dynare\5.5\matlab');
dynare soe_pac_full.mod
save('my_soe_pac_full.mat');

