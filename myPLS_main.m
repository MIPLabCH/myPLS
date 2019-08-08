% This is the main script of to run the myPLS Toolbox
% 
% ------------------------------STEPS-------------------------------------
%
% The script includes: 
%   1. Call of a script with PLS inputs and their description
%   2. Call of the functions to run PLS and plot the results
%
%
% --------------------------REQUIREMENTS----------------------------------
%
% Requires:
%   - SPM for saving results onto volume
%   - Slover to display slice maps
%   - Function ploterr (Copyright (c) 2008, Felix Zoergiebel) for bar plots
%
% ------------------------- ON PLS METHOD---------------------------------
% For general descriptions of PLS for medical image analysis, we refer to:
%
% ﻿- Krishnan, A., Williams, L.J., McIntosh, A.R., Abdi, H., 2011. Partial 
%    Least Squares (PLS) methods for neuroimaging: A tutorial and review. 
%    Neuroimage 56, 455–475. doi:10.1016/j.neuroimage.2010.07.034
%
%  - McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of 
%    neuroimaging data: Applications and advances. Neuroimage 23, 250–263. 
%    doi:10.1016/j.neuroimage.2004.07.020
%
%
% ----------------------------EXAMPLES------------------------------------
% Example applications, which used this toolbox:
%
%  - Behavior PLS with brain network-based measures:
%
%     * Zöller, D., Sandini, C., Karahanoğlu, F.I., Padula, M.C., Schaer, 
%       M., Eliez, S., Van De Ville, D., 2019. Large-scale brain network 
%       dynamics provide a measure of psychosis and anxiety in 22q11.2 
%       deletion syndrome. Biol. Psychiatry Cogn. Neurosci. Neuroimaging 
%       in press. doi:10.1016/j.bpsc.2019.04.004
%
%  - Behavior PLS with connectivity measures and different grouping in PLS
%    and resampling: 
% 
%     * Kebets, V., Holmes, A.J., Orban, C., Tang, S., Li, J., Sun, N.,
%       Kong, R., Poldrack, R.A., Yeo, B.T.T., 2019. Somatosensory-Motor
%       Dysconnectivity Spans Multiple Transdiagnostic Dimensions of
%       Psychopathology. Biol. Psychiatry.
%       doi:10.1016/j.biopsych.2019.06.013
%
%  - Contrast PLS for multivariate analysis of group differences and
%    developmental effects: 
%
%     * Zöller, D., Schaer, M., Scariati, E., Padula, M.C., Eliez, S., Van 
%       De Ville, D., 2017. Disentangling resting-state BOLD variability 
%       and PCC functional connectivity in 22q11.2 deletion syndrome. 
%       Neuroimage 149, 85–97. doi:10.1016/j.neuroimage.2017.01.064
%
%     * Zöller, D., Padula, M.C., Sandini, C., Schneider, M., Scariati, 
%       E., Van De Ville, D., Schaer, M., Eliez, S., 2018. Psychotic 
%       symptoms influence the development of anterior cingulate BOLD 
%       variability in 22q11.2 deletion syndrome. Schizophr. Res. 193, 
%       319–328. doi:10.1016/j.schres.2017.08.003
%
%
% -----------------------------CREDITS------------------------------------
% Code written by Prof. Dimitri Van De Ville, Daniela Zoeller and Valeria
% Kebets, with subfunctions borrowed from PLS toolbox by Rotman Baycrest
% (https://www.rotman-baycrest.on.ca/index.php?section=84)
% 
% These scripts and functions are based on myPLS scripts previously 
% published at https://miplab.epfl.ch/index.php/software/PLS
%
% 
% Please cite the following papers when using this code:
% 
% Zöller, D., Schaer, M., Scariati, E., Padula, M.C., Eliez, S., Van De 
% Ville, D., 2017. Disentangling resting-state BOLD variability and PCC 
% functional connectivity in 22q11.2 deletion syndrome. Neuroimage 149, 
% 85–97. doi:10.1016/j.neuroimage.2017.01.064
%
% Kebets, V., Holmes, A.J., Orban, C., Tang, S., Li, J., Sun, N., Kong, R.,
% Poldrack, R.A., Yeo, B.T.T., 2019. Somatosensory-Motor Dysconnectivity
% Spans Multiple Transdiagnostic Dimensions of Psychopathology. Biol.
% Psychiatry. doi:10.1016/j.biopsych.2019.06.013
% 
% McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of 
% neuroimaging data: Applications and advances. Neuroimage 23, 250–263. 
% doi:10.1016/j.neuroimage.2004.07.020

clc; clear all; close all;

addpath('./myPLS_functions')
addpath('./RotmanBaycrest')
addpath('./misc')
addpath('./spm12')


%% define all the inputs
% Modify this script to setup your PLS analysis
myPLS_inputs


%% check all inputs for validity
% !!! always run this function to check your setup before running PLS !!!
[input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);

%% save & plot input data
myPLS_plot_inputs(input,pls_opts,save_opts)

%% run PLS analysis, including permutation testing and bootstrapping
res = myPLS_analysis(input,pls_opts);

%% save & plot results data
% If you run multiple PLS analyses, correct the resulting p-values for
% multiple comparisons before executing the following function
myPLS_plot_results(res,save_opts);


