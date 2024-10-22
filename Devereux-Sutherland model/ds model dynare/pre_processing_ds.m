% Simple script to do the pre-processing of the model
% With Dynare, you need to execute this code once at the start of each
% session.
%Also, add the folder "common functions" to the MATLAB path.
%
% Copyright (C) 2024 Guillermo Hausmann Guil

clear;


addpath('C:\dynare\5.2\matlab');
dynare prepare_model_ds.mod
save('my_ds_model.mat');

