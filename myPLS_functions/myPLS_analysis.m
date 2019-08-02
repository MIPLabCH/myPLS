function [U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_img_loadings] = myPLS_analysis(X,Y,normalization_img,normalization_behav,nGroups_PLS,grouping_PLS)
%
% PLS analysis (main script)
%
% Inputs:
% - X                   : N x M matrix, N is #subjects, M is #imaging variables
% - Y                   : N x B matrix, B is #behaviors
% - normalization_img   : normalization options for imaging data
% - normalization_behav : normalization options for behavior data
%                         0 = no normalization
%                         1 = zscore across all subjects
%                         2 = zscore within groups (default)
%                         3 = std normalization across subjects (no centering)
%                         4 = std normalization within groups (no centering)
% - nGroups_PLS          : number of groups for PLS analysis, will determine
%                          if each group has its own set of saliences
%                          1 = the cross-covariance matrix is computed across all subjects
%                          2 = the cross-covariance matrix is computed within each group, and
%                          each group has its set of behavior saliences 
% - grouping_PLS         : N x 1 vector, subject grouping for PLS analysis
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
%
% Outputs:
% - U                    : B x L matrix, L is #latent components (LC), behavior saliences
% - S                    : L x L matrix, singular values
% - V                    : M x L matrix, imaging saliences
% - Lx                   : N x L matrix, imaging scores
% - Ly                   : N x L matrix, behavior scores
% - explCovLC            : covariance explained by each LC
% - LC_behav_loadings    : B x L matrix, behavior loadings
% - LC_img_loadings      : M x L matrix, imaging loadings
%

% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
nBehav = size(Y,2);

% Data normalization
X = myPLS_norm(X,nGroups_PLS,grouping_PLS,normalization_img);
Y = myPLS_norm(Y,nGroups_PLS,grouping_PLS,normalization_behav);

% Cross-covariance matrix
R = myPLS_cov(X,Y,nGroups_PLS,grouping_PLS);

% Singular value decomposition
[U,S,V] = svd(R,'econ');
nLCs = min(size(S)); % Number of latent components (LC)

% ICA convention: turn LCs such that max is positive
for iter_lc = 1:nLCs
    [~,idx] = max(abs(V(:,iter_lc)));
    if sign(V(idx,iter_lc))<0
        V(:,iter_lc) = -V(:,iter_lc);
        U(:,iter_lc) = -U(:,iter_lc);
    end
end

% Amount of covariance explained by each LC
explCovLC = (diag(S).^2) / sum(diag(S.^2));

% Imaging & behavioral composite scores
Lx = X * V;

if nGroups_PLS == 1
    
    Ly = Y * U;
    
elseif nGroups_PLS == 2
    Ly = nan(size(Lx));
    
    iter = 1;
    for iter_group = 1:nGroups_PLS        
        Usel = U(iter:iter + nBehav - 1,:);        
        for iter_group2 = 1:nGroups_PLS
            idx = find(grouping_PLS == iter_group2);
            Ysel = Y(idx,:);
            Lyy(iter_group,idx,:) = Ysel * Usel;
        end        
        iter = iter + nBehav;
    end
    
    for iter_group = 1:nGroups_PLS
        idx = find(grouping_PLS == iter_group);
        first = idx(1);
        last = idx(end);
        Ly(first:last,:) = Lyy(iter_group,first:last,:);
    end
    
    clear Lyy Usel Ysel idx iter first last
end

% Imaging loadings (Pearson's correlations between Lx and X)
for iter_lc = 1:nLCs
    for iter_img = 1:size(X,2)
        tmpy = Lx(:,iter_lc);
        tmpx = X(:,iter_img);
        r = corrcoef(tmpx,tmpy.');
        LC_img_loadings(iter_img,iter_lc) = r(1,2);
        clear tmpy tmpx r
    end
end

% Behavior loadings (Pearson's correlations between Ly and Y)
for iter_lc = 1:nLCs
    for iter_behav = 1:nBehav
        tmpy = Ly(:,iter_lc);
        tmpx = Y(:,iter_behav);
        r = corrcoef(tmpx,tmpy.');
        LC_behav_loadings(iter_behav,iter_lc) = r(1,2);
        clear tmpy tmpx r
    end
end