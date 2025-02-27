% Simple script to do the pre-processing of the growth model.
% With Dynare, you need to execute this code once at the start of each
% session.
% Also, add the folder "common functions" to the MATLAB path.

clear;
addpath('C:\dynare\5.5\matlab');
dynare prepare_model_growth.mod
save('my_growth_model.mat');

