% Script to define all Inputs for PLS for Medical Image Processing Toolbox
% 
% Parameters to be defined in this script are:
%
%   - input : struct containing input data for the analysis
%       - .X             : N x M matrix, N is #subjects, M is #imaging variables
%       - .Y             : N x B matrix, B is #behaviors
%       - .grouping  : N x 1 vector, subject grouping for PLS analysis
%                               e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                               subject 3 belongs to group 2.
%       - [.group_names]: Names of the groups (optional)
%       - [.behav_names]: Names of the behavior variables (optional) 
%   - pls_opts : options for the PLS analysis
%       - .nPerms              : number of permutations to run
%       - .nBootstraps         : number of bootstrapping samples to run
%       - .normalization_img   : normalization options for imaging data
%       - .normalization_behav : normalization options for behavior data
%              0 = no normalization
%              1 = zscore across all subjects
%              2 = zscore within groups (default)
%              3 = std normalization across subjects (no centering)
%              4 = std normalization within groups (no centering)
%       - .grouped_PLS         : binary variable indicating if groups
%                                should be considered when computing R
%              0 = PLS will computed over all subjects
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices ( as in conventional behavior PLS)
%       - .grouped_perm : binary variable indicating if groups should be 
%               considered during the permutations
%              0 = permutations ignoring grouping
%              1 = permutations within group
%       - .grouped_boot : binary variable indicating if groups should be 
%               considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - .boot_procrustes_mod : mode for bootstrapping procrustes transform
%              1 = standard (rotation computed only on U)
%              2 = average rotation of U and V
%       - [.behav_type]        : Type of behavioral analysis
%              'behavior' for standard behavior PLS
%              'contrast' to simply compute contrast between two groups
%              'contrastBehav' to combine contrast and behavioral measures)
%              'contrastBehavInteract' to also consider group-by-behavior interaction effects
%   - save_opts: Options for result saving and plotting
%       - .output_path   : path where to save the results
%       - [.prefix]      : prefix of all results files (optional)
%       - .img_type      : Specify how to plot the results
%              'volume' for voxel-based data in nifti Format - results 
%                       will be displayed as bootstrap ratios in a brain map
%              'corrMat' for ROI-to-ROI correlation matrix - results will 
%                       be displayed as bootstrap ratios in a correlation matrix
%              'barPlot' for any type of brain data in already vectorized 
%                       form - results will be displayed as barplots
%       - .mask_file     : gray matter mask, only required if imagingType='volume'
%       - .grouped_plots : binary variable indicating if groups should be 
%                          considered during plotting
%              0 = plotting ignoring grouping
%              1 = plotting cosidering grouping
%       - [.alpha]         : significance level for LCs [default = 0.05]




%% ---------- Input data ----------

%%%%%%%% LOAD YOUR DATA HERE %%%%%%%%
load('example_data.mat');
X0=brain_data;
Y0behav=[age,sex,FSIQ]; % sex - 0=female/1=male


% if you want to regress nuisance variables, do it here


% --- brain data ---
% Matrix X0 is typically a matrix with brain imaging data,
%  of size subjects (rows) x imaging features (columns)
input.brain_data=X0;

% --- behavior data ---
% Matrix behavData is a a matrix containing behavior data,
%  of size subjects (rows) x behavior features (columns)
%  Y0 will be constructed depending on the pls_opts.behav_type:
%   * if behav_type = 'contrast': behavData can be empty, Y0 will only depend
%      on the grouping information
%   * if behav_type = 'behavior': Y0 will be identical to behavData 
%   * if behav_type = 'contrastBehav', or 'contrastBehavInteract': 
%      Y0 will be constructed based on the grouping information and behavData
input.behav_data=Y0behav;

% --- grouping data ---
% subj_grouping: group assignment vector
% binary variable indicating the group, can contain multiple groups
input.grouping=diagnosis;


% --- Names of the groups ---
% here you can specify the names of the groups for the plots
% input.group_names={'group 1','group 2'}; 

% --- Names of the behavior data ---
input.behav_names={'age','sex','IQ'}; 



clear age diagnosis sex FSIQ brain_data

%% ---------- Options for PLS ----------

% --- Permutations & Bootstrapping ---
pls_opts.nPerms = 200;
pls_opts.nBootstraps = 100;

% --- Data normalization options ---
% 0: no normalization
% 1: zscore across all subjects (default)
% 2: zscore within groups
% 3: std normalization across subjects (no centering)
% 4: std normalization within groups (no centering)
pls_opts.normalization_img=1;
pls_opts.normalization_behav=1;

% --- PLS grouping option ---
% 0: PLS will computed over all subjects
% 1: R will be constructed by concatenating group-wise covariance matrices
%     (as in conventional behavior PLS, see Krishnan et al., 2011)
pls_opts.grouped_PLS=1; 


% --- Permutations grouping option ---
% 0: permutations ignoring grouping
% 1: permutations within group
pls_opts.grouped_perm=0;

% --- Bootstrapping grouping option ---
% 0: bootstrapping ignoring grouping
% 1: bootstrapping within group
pls_opts.grouped_boot=0;

% --- Mode for bootstrapping procrustes transform ---
% in some cases, rotation only depending on U results in extremely low
% standard errors and bootstrap ratios close to infinity
% in mode 2, we therefore compute the transformation matrix both on U and V
pls_opts.boot_procrustes_mod=1; % 1 - standard; 2 - average rotation of U and V

% --- Type of behavioral analysis ---
% 'behavior' for standard behavior PLS
% 'contrast' to simply compute contrast between two groups
% 'contrastBehav' to combine contrast and behavioral measures
% 'contrastBehavInteract' to also consider group-by-behavior interaction effects
pls_opts.behav_type = 'behavior';


%% ---------- Options for result saving and plotting ----------
% --- path where to save the results ---
save_opts.output_path='./example_results';

% --- prefix of all results files ---
% this is also the default prefix of the toolbox if you don't define
% anything
save_opts.prefix = sprintf('myPLS_TYPE%s_NORM%d%d',pls_opts.behav_type,...
    pls_opts.normalization_img, pls_opts.normalization_behav);


% --- Type of brain data ---
% Specify how to plot the results
% 'volume'  for voxel-based data in nifti Format - results will be 
%           displayed as bootstrap ratios in a brain map
% 'corrMat' for ROI-to-ROI correlation matrix - results will be displayed
%           as bootstrap ratios in a correlation matrix
% 'barPlot' for any type of brain data in already vectorized form - results 
%           will be displayed as barplots
save_opts.img_type = 'volume';

% --- Brain mask ---
% (gray matter mask, only required if imagingType='volume')
save_opts.mask_file = 'example_mask.nii'; % filename of binary mask that will constrain analysis


% --- Plotting grouping option ---
% 0: Plots ignoring grouping
% 1: Plots considering grouping
save_opts.grouped_plots=1;


% --- Significance level for latent components ---
save_opts.alpha=0.1; % for the sake of the example data




