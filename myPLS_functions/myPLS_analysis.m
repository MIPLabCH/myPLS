function res = myPLS_analysis(input,pls_opts)
%
% PLS analysis (main script)
%
% Inputs:
%   - input : struct containing input data for the analysis
%       - .brain_data    : N x M matrix, N is #subjects, M is #imaging variables
%       - .behav_data    : N x B matrix, B is #behaviors
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
%       - [.grouped_PLS]       : binary variable indicating if groups
%                                should be considered when computing R
%              0 = PLS will computed over all subjects [only option for contrast PLS] 
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices [default for behavior PLS]
%       - [.grouped_perm] : binary variable indicating if groups should be 
%               considered during the permutations
%              0 = permutations ignoring grouping
%              1 = permutations within group
%       - [.grouped_boot] : binary variable indicating if groups should be 
%               considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - [.boot_procrustes_mod]: mode for bootstrapping procrustes transform
%              1 = standard (rotation computed only on U) [default]
%              2 = average rotation of U and V
%       - [.behav_type]        : Type of behavioral analysis
%              'behavior' for standard behavior PLS [default]
%              'contrast' to simply compute contrast between two groups
%              'contrastBehav' to combine contrast and behavioral measures)
%              'contrastBehavInteract' to also consider group-by-behavior interaction effects

%
% Outputs:
%   res: Struct containing all results
%     .X0, .Y0     : unnormalized input matrices
%     .X, .Y       : normalized input matrices
%     .U           : B x L matrix, L is #latent components (LC), behavior saliences
%     .V           : M x L matrix, imaging saliences
%     .S           : L x L matrix, singular values (diagonal matrix)
%     .explCovLC   : covariance explained by each LC
%     .myLVpvals   : p-value for each LV (from permutation testing
%     .Ub_vect     : 3D matrix with bootstrapping samples of U
%     .Vb_vect     : 3D matrix with bootstrapping samples of V
%


%% initialize
% this function checks all relevant inputs for validity
[input,pls_opts] = myPLS_analysis_initialize(input,pls_opts);


% number of subjects
nSubj = size(input.brain_data,1); 

% number of behavior scores
nBehav=size(input.behav_data,2);

% number of imaging measures
nImg = size(input.brain_data,2);  

% number and IDs of groups
groupIDs=unique(input.grouping);


%% get brain matrix X0
X0=input.brain_data;


%% construct the behavior/contrast/interaction matrix Y0
[Y0,design_names,nDesignScores] = myPLS_getY(pls_opts.behav_type,...
    input.behav_data,input.grouping,input.group_names,input.behav_names);

disp('... Input data information ...')
disp(['Number of observations (subjects): ' num2str(nSubj)]);
disp(['Number of brain measures (voxels/connections): ' num2str(nImg)]);
disp(['Number of behavior measures: ' num2str(nBehav)]);
disp(['Number of design measures (behavior/contrasts): ' num2str(nDesignScores)]);
disp(' ')


%% Normalize input data X and Y 
% (if there is a group contrast in Y, groups won't be considered in any case)
X = myPLS_norm(X0,input.grouping,pls_opts.normalization_img);
Y = myPLS_norm(Y0,input.grouping,pls_opts.normalization_behav);

%% Cross-covariance matrix
R = myPLS_cov(X,Y,input.grouping,pls_opts.grouped_PLS);

%% Singular value decomposition
[U,S,V] = svd(R,'econ');
nLCs = min(size(S)); % Number of latent components (LC)

% ICA convention: turn LCs such that max is positive
for iLC = 1:nLCs
    [~,maxID] = max(abs(V(:,iLC)));
    if sign(V(maxID,iLC))<0
        V(:,iLC) = -V(:,iLC);
        U(:,iLC) = -U(:,iLC);
    end
end

% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

% Compute Imaging & behavioral composite scores
[Lx,Ly,LC_img_loadings,LC_behav_loadings] = myPLS_get_PLSscores(X,Y,V,U,input.grouping,pls_opts);

%% Permutation testing for LV significance
% !!! permutations should be run using the already normalized X and Y !!!
Sp_vect = myPLS_permutations(X,Y,U,input.grouping,pls_opts);

% compute the p-values from the permutation null distribution
LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts);

%% Bootstrapping to test stability of brain saliences
% !!! use non-normalized X0 and Y0, normalization will be done again because of resampling WITH replacement !!!
boot_results=myPLS_bootstrapping(X0,Y0,U,V,input.grouping,pls_opts);


%% save all result variables in struct
res.X0=X0;
res.Y0=Y0;
res.X=X;
res.Y=Y;
res.design_names=design_names;
res.R=R;
res.U=U;
res.S=S;
res.V=V;
res.explCovLC=explCovLC;
res.Lx=Lx;
res.Ly=Ly;
res.LC_img_loadings=LC_img_loadings;
res.LC_behav_loadings=LC_behav_loadings;
%res.Sp_vect=Sp_vect; % do we need to save this?
res.LC_pvals=LC_pvals;
res.boot_results=boot_results;













