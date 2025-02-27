% Simple script to do the pre-processing of the SOE model (calibration for
% comparative statics).
% With Dynare, you need to execute this code once at the start of each
% session.
% Also, add the folder "common functions" to the MATLAB path.

clear;
addpath('C:\dynare\5.5\matlab');
dynare soe_uzawa_basic.mod
save('my_soe_uzawa_basic.mat');

