%% This is the main script to run myPLS Toolbox
% 
% ------------------------------STEPS-------------------------------------
%
% The script includes: 
%   1. Call of a script with PLS inputs and their description
%   2. Call of the functions to run PLS and plot the results


clc; clear all; close all;

addpath('./myPLS_functions')
addpath('./RotmanBaycrest')
addpath('./misc')
%addpath('./spm12')

%% Define all the inputs
% Modify this script to setup your PLS analysis

myPLS_inputs

%% Check all inputs for validity
% !!! always run this function to check your setup before running PLS !!!

[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);

%% Save & plot input data

% myPLS_plot_inputs(input,pls_opts,save_opts)

%% Run PLS analysis (including permutation testing and bootstrapping)

res = myPLS_analysis(input,pls_opts);

%% Save & plot results
% If you run multiple PLS analyses, correct the resulting p-values for
% multiple comparisons before executing the following function

myPLS_plot_results(res,save_opts);
