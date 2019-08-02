% Script to define all Inputs for PLS for Medical Image Processing Toolbox
% 
% Parameters to be defined in this script are:
%
%   input: struct containing input data for the analysis
%     .brainData: brain data (nSub x nData1) - e.g. voxelwise seed connectivity
%     .behavData: behavior data (nSub x nData2) - e.g. [age, IQ]
%     .subj_grouping: binary grouping information (nSub x 1)
%   plsOpts: Options for the PLS analysis
%     .behavType: Type of behavioral analysis
%        'behavior' for standard behavior PLS
%        'contrast' to simply compute contrast between two groups
%        'contrastBehav' to combine contrast and behavioral measures)
%        'contrastBehavInteract' to also consider group-by-behavior interaction effects
%     .NUM_PERMS: number of permutations to run
%     .NUM_BOOTSTRAP: number of bootstrapping samples to run
%     .CONST_NORM_IMAGING and .CONST_NORM_BEHAV: Data normalization options
%        0: no normalization
%        1: zscore across all subjects (default) 
%        2: zscore within groups
%        3: std normalization across subjects (no centering)
%        4: std normalization within groups (no centering)
%   saveOpts: Options for result saving and plotting
%     .CONST_OUTPUT_PATH: path where to save the results
%     [.prefix]: prefix of all results files (optional)
%     [.CONST_GROUPS]: Names of the groups (optional)
%     [.CONST_BEHAV_NAMES]: Names of the behavior variables (optional) 
%     .imagingType: Specify how to plot the results
%        'volume'  for voxel-based data in nifti Format - results will be 
%           displayed as bootstrap ratios in a brain map
%        'corrMat' for ROI-to-ROI correlation matrix - results will be 
%           displayed as bootstrap ratios in a correlation matrix
%        'barPlot' for any type of brain data in already vectorized form 
%           - results will be displayed as barplots
%     .maskFile - gray matter mask, only required if imagingType=volume




%% ---------- Input data ----------

%%%%%%%% LOAD YOUR DATA HERE %%%%%%%%
load('example_data.mat');
X0=brain_data;
Y0behav=[age,sex,FSIQ]; % sex - 0=female/1=male


% if you want to regress nuisance variables, do it here


% --- brain data ---
% Matrix X0 is typically a matrix with brain imaging data,
%  of size subjects (rows) x imaging features (columns)
input.brainData=X0;

% --- behavior data ---
% Matrix Y0behav is a a matrix containing behavior data,
%  of size subjects (rows) x behavior features (columns)
%  Y0 will be constructed depending on the plsOpts.behavType:
%   * if behavType = 'contrast': Y0behav can be empty, Y0 will only depend
%      on the grouping information
%   * if behavType = 'behavior': Y0 will be identical to Y0behav 
%   * if behavType = 'contrastBehav', or 'contrastBehavInteract': 
%      Y0 will be constructed based on the grouping information and Y0behav
input.behavData=Y0behav;

% --- grouping data ---
% subj_grouping: group assignment vector
% binary variable indicating the group, can contain multiple groups
input.subj_grouping=diagnosis;


clear age diagnosis sex FSIQ brain_data

%% ---------- Options for PLS ----------

% --- Type of behavioral analysis ---
% 'behavior' for standard behavior PLS
% 'contrast' to simply compute contrast between two groups
% 'contrastBehav' to combine contrast and behavioral measures
% 'contrastBehavInteract' to also consider group-by-behavior interaction effects
plsOpts.behavType = 'behavior';

% --- Permutations & Bootstrapping ---
plsOpts.NUM_PERMS = 200;
plsOpts.NUM_BOOTSTRAP = 100;

% --- Data normalization options ---
% 0: no normalization
% 1: zscore across all subjects (default)
% 2: zscore within groups
% 3: std normalization across subjects (no centering)
% 4: std normalization within groups (no centering)
plsOpts.CONST_NORM_IMAGING=1;
plsOpts.CONST_NORM_BEHAV=1;


% mode for bootstrapping procrustes transform
% in some cases, rotation only depending on U results in extremely low
% standard errors and bootstrap ratios close to infinity
% in mode 2, we therefore compute the transformation matrix both on U and V
plsOpts.PROCRUSTES_MOD=1; % 1 - standard; 2 - average rotation of U and V


%% ---------- Options for result saving and plotting ----------
% --- path where to save the results ---
saveOpts.CONST_OUTPUT_PATH='/name/of/path';

% --- prefix of all results files ---
saveOpts.prefix = 'myPLSresults_example';


% --- Names of the groups ---
% here you can specify the names of the groups for the plots
saveOpts.CONST_GROUPS={'group 1','group 2'}; 

% --- Names of the behavior data ---
saveOpts.CONST_BEHAV_NAMES={'age','sex','IQ'}; 



% --- Type of brain data ---
% Specify how to plot the results
% 'volume'  for voxel-based data in nifti Format - results will be 
%           displayed as bootstrap ratios in a brain map
% 'corrMat' for ROI-to-ROI correlation matrix - results will be displayed
%           as bootstrap ratios in a correlation matrix
% 'barPlot' for any type of brain data in already vectorized form - results 
%           will be displayed as barplots
saveOpts.imagingType = 'volume';

% --- Brain mask ---
% (gray matter mask, only required if imagingType='volume')
saveOpts.maskFile = fullfile('path','to','maskFile.nii'); % filename of binary mask that will constrain analysis 







