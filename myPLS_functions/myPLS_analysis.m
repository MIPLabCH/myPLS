function res = myPLS_analysis(input,pls_opts)

% This is the main function computing the PLS analysis
%
% Inputs:
% - input : struct containing input data for the analysis
%       - .brain_data  : N x M matrix, N is #subjects, M is #imaging variables
%       - .behav_data  : N x B matrix, B is #behavior/design variables
%       - .grouping    : N x 1 vector, subject grouping for PLS analysis
%                        e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                        subject 3 belongs to group 2.
%       - .group_names : names of the groups (optional)
%       - .behav_names : names of the behavior variables (optional) 
% - pls_opts : options for the PLS analysis
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
%              0 = PLS will computed over all subjects [only option for contrast PLS] 
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices [default for behavior PLS]
%       - .grouped_perm        : binary variable indicating if groups
%                                should be considered during permutations
%              0 = permutations ignoring grouping
%              1 = permutations within group
%       - .grouped_boot        : binary variable indicating if groups
%                                should be considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - .boot_procrustes_mod : mode for Procrustes transform during bootstrapping
%              1 = standard (rotation computed only on U) [default]
%              2 = average rotation of U and V
%       - .save_boot_resampling : binary variable indicating if bootstrap
%                                 resampling data should be saved or not
%              0 = no saving of bootstrapping resampling data
%              1 = save bootstrapping resampling data
%       - .behav_type           : type of behavior analysis
%              'behavior' for standard behavior PLS [default]
%              'contrast' to simply compute contrast between two groups
%              'contrastBehav' to combine contrast and behavioral measures)
%              'contrastBehavInteract' to also consider group-by-behavior interaction effects
%
% Outputs:
%   - res: struct containing all results
%     - .X0, .Y0          : not normalized input matrices
%     - .X, .Y            : normalized input matrices
%     - .design_names     : names of design variables (changed only for contrast PLS)
%     - .grouping,.group_names : copied from input for plotting
%     - .R                : L x M matrix, L is #latent components (LC),
%                           brain-behavior covariance matrix
%     - .U                : B x L matrix, behavior saliences 
%     - .V                : M x L matrix, imaging saliences
%     - .S                : L x L matrix, singular values (diagonal matrix)
%     - .explCovLC        : covariance explained by each LC
%     - .LC_pvals         : p-value for each LC (from permutation testing)
%     - .Lx               : brain scores
%     - .Ly               : behavior/design scores
%     - .LC_img_loadings  : corr(Lx,X)
%     - .LC_behav_loadings: corr(Ly,Y)
%     - .boot_results     : struct with bootstrapping results
%           - .Ub_vect    : 3D matrix with bootstrapping samples of U
%           - .Vb_vect    : 3D matrix with bootstrapping samples of V
%           - .Lxb,.Lyb,.LC_img_loadings_boot,.LC_behav_loadings_boot :
%                           3D matrices with bootstrapping samples of PLS scores
%           - .*_mean     : mean of bootstrapping distributions
%           - .*_std      : standard deviation of bootstrapping distributions
%           - .*_lB       : lower bound of 95% confidence interval of bootstrapping distributions
%           - .*_uB       : upper bound of 95% confidence interval of bootstrapping distributions

%% Get constants

% Number of subjects
nSubj = size(input.brain_data,1); 

% Number of behavior/design variables
nBehav = size(input.behav_data,2);

% Number of imaging variables
nImg = size(input.brain_data,2);  

%% Get imaging matrix X0

X0 = input.brain_data;

%% Generate the behavior/contrast/interaction matrix Y0

[Y0,design_names,nDesignScores] = myPLS_getY(pls_opts.behav_type,...
    input.behav_data,input.grouping,input.group_names,input.behav_names);

disp('... Input data information ...')
disp(['Number of observations (subjects): ' num2str(nSubj)]);
disp(['Number of brain measures (voxels/connections): ' num2str(nImg)]);
disp(['Number of design measures (behavior/contrasts): ' num2str(nDesignScores)]);
disp(' ')

%% Normalize input data matrices X and Y 
% (if there is a group contrast in Y, groups won't be considered in any case)

[X,meanX0,stdX0] = myPLS_norm(X0,input.grouping,pls_opts.normalization_img);
[Y,meanY0,stdY0] = myPLS_norm(Y0,input.grouping,pls_opts.normalization_behav);

%% Cross-covariance matrix
R = myPLS_cov(X,Y,input.grouping,pls_opts.grouped_PLS);

%% Singular value decomposition
[U,S,V] = svd(R,'econ');

% Number of latent components (LCs)
nLC = min(size(S)); 

% ICA convention: turn LCs such that max is positive
for iLC = 1:nLC
    [~,maxID] = max(abs(V(:,iLC)));
    if sign(V(maxID,iLC))<0
        V(:,iLC) = -V(:,iLC);
        U(:,iLC) = -U(:,iLC);
    end
end

% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

% Compute PLS scores & loadings
[Lx,Ly,corr_Lx_X,corr_Ly_Y,corr_Lx_Y,corr_Ly_X] = ...
    myPLS_get_PLS_scores_loadings(X,Y,V,U,input.grouping,pls_opts);

%% Permutation testing for LC significance

Sp_vect = myPLS_permutations(X,Y,U,input.grouping,pls_opts);

% Compute the p-values from the permutation null distribution
LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts);

%% Bootstrapping to test stability of PLS loadings

boot_results = myPLS_bootstrapping(X0,Y0,U,V,S,input.grouping,pls_opts);

%% Save all result variables in struct

res.X0 = X0;
res.meanX0=meanX0;
res.stdX0=stdX0;
res.Y0 = Y0;
res.meanY0=meanY0;
res.stdY0=stdY0;
res.X = X;
res.Y = Y;
res.design_names = design_names;
res.grouping = input.grouping;
res.group_names = input.group_names;
if isfield(input,'img_names'); res.img_names = input.img_names; end
res.R = R;
res.U = U;
res.S = S;
res.V = V;
res.explCovLC = explCovLC;
res.LC_pvals = LC_pvals;
res.Lx = Lx;
res.Ly = Ly;
res.LC_img_loadings = corr_Lx_X; %%% can be changed to corr_Ly_X
res.LC_behav_loadings = corr_Ly_Y; %%% can be changed to corr_Lx_Y
res.Sp_vect = Sp_vect;
res.LC_pvals = LC_pvals;
res.boot_results = boot_results;
