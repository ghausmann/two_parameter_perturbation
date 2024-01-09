% Simple script to do the pre-processing of the SOE model (calibration with
% Mexican data).
% With Dynare, you need to execute this code once at the start of each
% session.
%Also, add the folder "common functions" to the MATLAB path.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

clear;
addpath('C:\dynare\5.2\matlab');
dynare soe_uzawa_full.mod
save('my_soe_uzawa_full.mat');

