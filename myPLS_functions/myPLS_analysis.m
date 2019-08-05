function res = myPLS_analysis(input,pls_opts)
%
% PLS analysis (main script)
%
% Inputs:
%   - input : struct containing input data for the analysis
%       - .brain_data    : N x M matrix, N is #subjects, M is #imaging variables
%       - .behav_data    : N x B matrix, B is #behaviors
%       - .grouping_PLS  : N x 1 vector, subject grouping for PLS analysis
%                               e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                               subject 3 belongs to group 2.
%       - [.group_names]: Names of the groups (optional)
%       - [.behav_names]: Names of the behavior variables (optional) 
%   - pls_opts : options for the PLS analysis
%       - [.behav_type]          : Type of behavioral analysis
%              'behavior' for standard behavior PLS
%              'contrast' to simply compute contrast between two groups
%              'contrastBehav' to combine contrast and behavioral measures)
%              'contrastBehavInteract' to also consider group-by-behavior interaction effects
%       - .nPerms              : number of permutations to run
%       - .nBootstraps         : number of bootstrapping samples to run
%       - .normalization_img   : normalization options for imaging data
%       - .normalization_behav : normalization options for behavior data
%              0 = no normalization
%              1 = zscore across all subjects
%              2 = zscore within groups (default)
%              3 = std normalization across subjects (no centering)
%              4 = std normalization within groups (no centering)

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
% default: behavior PLS
if ~isfield(pls_opts,'behav_type') || isempty(pls_opts.behav_type)
    pls_opts.behav_type='behavior';
end

% compatibility with X0 and Y0 inputs
if ~isfield(input,'brain_data') && isfield(input,'X0')
    input.brain_data=input.X0;
elseif ~isfield(input,'brain_data') && ~isfield(input,'X0')
    error('brain data input missing');
end
if ~isfield(input,'behav_data') && isfield(input,'Y0')
    input.behav_data=input.Y0;
elseif ~isfield(input,'behav_data') && ~isfield(input,'Y0')
    error('behavior data input missing');
end 

% Check that dimensions of X & Y are correct
if(size(input.brain_data,1) ~= size(input.behav_data,1))
    error('Input arguments X and Y should have the same number of rows');
end

% number of subjects
nSubj = size(input.X0,1); 

% number of behavior scores
nBehav=size(input.behav_data,2);

% number of imaging measures
nImg = size(input.X0,2);  

% number and IDs of groups
groupIDs=unique(input.grouping_PLS);
nGroups=length(groupIDs);

% create defaults for group and behavior names, if not specified
if ~isfield(input,'group_names') || isempty(input.group_names)
    input.group_names=cell(nGroups,1);
    for iG=1:nGroups
        input.group_names{iG}=['group ' num2str(groupIDs(iG))];
    end
end
if ~isfield(input,'behav_names') || isempty(input.behav_names)
    input.behav_names=cell(nBehav,1);
    for iB=1:nBehav
        input.behav_names{iB}=['behavior ' num2str(iB)];
    end
end


%% get brain matrix X0
X0=input.brain_data;


%% construct the behavior/contrast/interaction matrix Y0
[Y0,design_names,nDesignScores] = myPLS_getY(pls_opts.behav_type,...
    input.behav_data,input.grouping_PLS,input.group_names,input.behav_names);

disp(['Number of observations (subjects): ' num2str(nSubj)]);
disp(['Number of brain measures (voxels/connections): ' num2str(nImg)]);
disp(['Number of behavior measures: ' num2str(nBehav)]);
disp(['Number of design measures (behavior/contrasts): ' num2str(nDesignScores)]);

%% set up grouping variable for contrast analyses
% TODO: here, we should replace the *_analysis variables by separate options
% to select whether grouping should be considered for PLS and/or
% bootstrapping separately

if contains(pls_opts.behav_type,'contrast')
    subj_grouping_analysis=ones(size(input.grouping_PLS));
    
    if pls_opts.normalization_behav==2 || pls_opts.normalization_behav==4 || ...
            pls_opts.normalization_img==2 || pls_opts.normalization_img==4
        warning('Normalization within groups selected, but contrast in Y -> normalization will be done across all subjects!')
    end
else
    subj_grouping_analysis=input.grouping_PLS;
end


%% Normalize input data X and Y 
% (if there is a group contrast in Y, groups won't be considered in any case)
X = myPLS_norm(input.X0,subj_grouping_analysis,pls_opts.normalization_img);
Y = myPLS_norm(Y0,subj_grouping_analysis,pls_opts.normalization_behav);

%% Cross-covariance matrix
R = myPLS_cov(X,Y,subj_grouping_analysis);

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


% %% Permutation testing for LV significance
% % !!! permutations should be run using the already normalized X and Y !!!
% Sp_vect = myPLS_permut(X,Y,U,subj_grouping_analysis,pls_opts);
% 
% % compute the p-values from the permutation null distribution
% myLCpvals = myPLS_getLCpvals(Sp_vect,S,pls_opts);
% 
% %% Bootstrapping to test stability of brain saliences
% % !!! use non-normalized X0 and Y0, normalization will be done again because of resampling WITH replacement !!!
% [Ub_vect,Vb_vect]=myPLS_bootstrapping(X0,Y0,U,subj_grouping_analysis,plsOpts);

%% save all result variables in struct
res.X0=X0;
res.Y0=Y0;
res.X=X;
res.Y=Y;
res.design_names=design_names;
res.nDesignScores=nDesignScores;
res.R=R;
res.U=U;
res.S=S;
res.V=V;
res.explCovLC=explCovLC;
% res.myLCpvals=myLCpvals;
% res.Ub_vect=Ub_vect;
% res.Vb_vect=Vb_vect;






% 
% % Imaging & behavioral composite scores
% Lx = X * V;
% 
% if nGroups_PLS == 1
%     
%     Ly = Y * U;
%     
% elseif nGroups_PLS == 2
%     Ly = nan(size(Lx));
%     
%     iter = 1;
%     for iter_group = 1:nGroups_PLS        
%         Usel = U(iter:iter + nBehav - 1,:);        
%         for iter_group2 = 1:nGroups_PLS
%             maxID = find(grouping_PLS == iter_group2);
%             Ysel = Y(maxID,:);
%             Lyy(iter_group,maxID,:) = Ysel * Usel;
%         end        
%         iter = iter + nBehav;
%     end
%     
%     for iter_group = 1:nGroups_PLS
%         maxID = find(grouping_PLS == iter_group);
%         first = maxID(1);
%         last = maxID(end);
%         Ly(first:last,:) = Lyy(iter_group,first:last,:);
%     end
%     
%     clear Lyy Usel Ysel idx iter first last
% end
% 
% % Imaging loadings (Pearson's correlations between Lx and X)
% for iLC = 1:nLCs
%     for iter_img = 1:size(X,2)
%         tmpy = Lx(:,iLC);
%         tmpx = X(:,iter_img);
%         r = corrcoef(tmpx,tmpy.');
%         LC_img_loadings(iter_img,iLC) = r(1,2);
%         clear tmpy tmpx r
%     end
% end
% 
% % Behavior loadings (Pearson's correlations between Ly and Y)
% for iLC = 1:nLCs
%     for iter_behav = 1:nBehav
%         tmpy = Ly(:,iLC);
%         tmpx = Y(:,iter_behav);
%         r = corrcoef(tmpx,tmpy.');
%         LC_behav_loadings(iter_behav,iLC) = r(1,2);
%         clear tmpy tmpx r
%     end
% end