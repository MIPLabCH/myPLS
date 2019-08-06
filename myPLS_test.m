% Testing script to check some of the PLS Toolbox functionalities
clc; clear all; close all;


addpath('./myPLS_functions')
addpath('./RotmanBaycrest')
addpath('./misc')


%% defining all the inputs
myPLS_inputs

%% run PLS analysis, including permutation testing and bootstrapping
res = myPLS_analysis(input,pls_opts);

% test setting up of default group & behavior names
input=rmfield(input,'group_names');
input=rmfield(input,'behav_names');
pls_opts=rmfield(pls_opts,'grouped_PLS');

% test comptibiltiy with X0 and Y0 input
input.X0=input.brain_data;
input.Y0=input.behav_data;
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 1 passed ----');
catch
    disp('---- Test 1 failed ----');
end


myPLS_inputs
% test of grouped PLS option
pls_opts.grouped_PLS=1;
try
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 2 passed ----');
catch
    disp('---- Test 2 failed ----');
end

% test grouping and contrast conflict
pls_opts.grouped_PLS=1;
pls_opts.behav_type='contrast';
try
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 3 failed: no error if conflicting grouping and contrast setup');
catch
    disp('---- Test 3 passed: error if conflicting grouping and contrast setup');
end

% test contrast PLS
pls_opts.grouped_PLS=0;
pls_opts.behav_type='contrast';
try
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 4 passed ----');
catch
    disp('---- Test 4 failed ----');
end

% test contrast & behavior PLS
pls_opts.grouped_PLS=0;
pls_opts.behav_type='contrastBehav';
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 5 passed ----');
catch
    disp('---- Test 5 failed ----');
end

pls_opts.grouped_PLS=0;
pls_opts.behav_type='contrastBehavInteract';
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 6 passed ----');
catch
    disp('---- Test 6 failed ----');
end

pls_opts.grouped_PLS=0;
pls_opts.behav_type='contrastBehavInter';
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 7 failed: no error if invalid behav_type');
catch
    disp('---- Test 7 passed: error if invalid behav_type');
end

pls_opts.behav_type='contrastBehavInteract';
pls_opts.boot_procrustes_mod=3;
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 8 failed: no error if invalid boot_procrustes_mod');
catch
    disp('---- Test 8 passed: error if invalid boot_procrustes_mod');
end

pls_opts.boot_procrustes_mod=2;
try 
    res = myPLS_analysis(input,pls_opts);
    disp('---- Test 9 passed ----');
catch
    disp('---- Test 9 failed ----');
end

