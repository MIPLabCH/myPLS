function [Lx,Ly,corr_Lx_X,corr_Ly_Y,corr_Lx_Y,corr_Ly_X] = ...
    myPLS_get_PLS_scores_loadings(X,Y,V,U,grouping,pls_opts)

% This function computes PLS scores & loadings.
% Individual-specific PLS scores are subjects' expression of imaging and 
% behavior/design saliences (covariance patterns).
% Component-specific PLS loadings express the contribution of original 
% imaging or behavior/design variables to the latent components (LCs).
%
% Inputs:
% - X           : N x M, N is #subjects, M is #imaging variables, imaging data
% - Y           : N x B, B is #behavior/design variables, behavior/design data
% - V           : M x L, L is #LCs, imaging saliences
% - U           : B x (L x #groups), behavior/design saliences
% - grouping    : N x 1 vector, subject grouping (e.g. diagnosis)
%                 e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                 subject 3 belongs to group 2
% - pls_opts    : options for the PLS analysis
% 
% Outputs:
% - Lx          : N x L, PLS imaging scores
% - Ly          : N x (L x #groups), PLS behavior/design scores
% - corr_Lx_X   : M x L, Pearson's correlation between imaging data and 
%                 PLS imaging scores (i.e., structure coefficients)
% - corr_Ly_Y   : B x (L x #groups), Pearson's correlation between behavior/design data
%                 and PLS behavior/design scores (i.e. structure coefficients)
% - corr_Lx_Y   : M x L, Pearson's correlation between 
%                 behavior/design data and PLS imaging scores
% - corr_Ly_X   : B x (L x #groups), Pearson's correlation between 
%                 behavior/design data and PLS imaging scores


if ~pls_opts.grouped_PLS
    grouping = ones(size(grouping));
end

% Number and IDs of groups
groupIDs = unique(grouping);
nGroups = length(groupIDs);

% Number of behavior variables
nBehav = size(Y,2);

% Number of imaging variables
%nImg = size(X,2);  

% Number of latent components (LCs)
nLC = size(V,2); 

%% Compute PLS scores

% PLS imaging  scores
Lx = X * V;

% PLS behavior/design scores
iter = 1;
for iG = 1:nGroups
    % Select U for this group
    idx = iter:iter + nBehav - 1;
    Usel = U(idx,:); 
    
    % Select Y for this group
    for iG2 = 1:nGroups
        this_groupID = find(grouping == groupIDs(iG2));
        Ysel = Y(this_groupID,:);
        Lyy(iG,this_groupID,:) = Ysel * Usel;
    end        
    iter = iter + nBehav;
end

% Stack group-specific behavior/design scores
for iG = 1:nGroups
    this_groupID = find(grouping == groupIDs(iG));
    if nLC==1
        Ly(this_groupID,:) = Lyy(iG,this_groupID,:)';
    else
        Ly(this_groupID,:) = Lyy(iG,this_groupID,:);
    end
end

%% Compute PLS loadings
 
% Correlations between Lx and X (imaging structure coefficients) 
corr_Lx_X = corr(Lx,X)';

% Correlations between Ly and Y (behavior/design structure coefficients) 
iter = 1;
for iG = 1:nGroups
    idx = iter:iter + nBehav - 1;
    this_groupID = find(grouping == groupIDs(iG));
    corr_Ly_Y(idx,:) = corr(Ly(this_groupID,:),Y(this_groupID,:))';
    iter = iter + nBehav;
end

% Correlations between Lx and Y
iter = 1;
for iG = 1:nGroups
    idx = iter:iter + nBehav - 1;
    this_groupID = find(grouping == groupIDs(iG));
    corr_Lx_Y(idx,:) = corr(Lx(this_groupID,:),Y(this_groupID,:))';
    iter = iter + nBehav;
end

% Correlations between Ly and X
corr_Ly_X = corr(Ly,X)';
