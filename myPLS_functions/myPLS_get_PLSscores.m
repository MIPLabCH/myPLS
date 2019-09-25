function [Lx,Ly,LC_img_loadings,LC_behav_loadings] = myPLS_get_PLSscores(X,Y,V,U,grouping,pls_opts)

% This function computes PLS scores & loadings
%
% Inputs:
% - X           : N x M, N is #subjects, M is #imaging variables, imaging data
% - Y           : N x B, B is #behavior/design variables, behavior/design data
% - V           : M x L, L is #latent components (LCs), imaging saliences
% - U           : B x (L x nGroups), behavior/design saliences
% - grouping    : N x 1 vector, subject grouping (e.g. diagnosis)
%                 e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                 subject 3 belongs to group 2
% - pls_opts    : options for the PLS analysis
% 
% Outputs:
% - Lx          : N x L, PLS imaging scores
% - Ly          : N x (L x nGroups), PLS behavior/design scores
% - corr_Lx_X   : M x L, Pearson's correlation between imaging data and 
%                 PLS imaging scores (structure coefficients)
% - corr_Ly_Y   : B x (L x nGroups), Pearson's correlation between behavior/design data
%                 and PLS behavior/design scores (structure coefficients)
% - corr_Lx_Y   : M x L, Pearson's correlation between 
%                 behavior/design data and PLS imaging scores
% - corr_Ly_X   : B x (L x nGroups), Pearson's correlation between 
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
nImg = size(X,1);  

% Number of latent components (LCs)
nLC = size(V,2); 

%% Compute PLS scores

% PLS imaging  scores
Lx = X * V;

% PLS behavior/design scores
% Dani: I removed the selection depending on the number of groups, I believe
% it should work like this even for one group
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
    Ly(this_groupID,:) = Lyy(iG,this_groupID,:);
end

%% Compute PLS loadings
 
% Pearson's correlations between Lx and X (imaging structure coefficients) 
% for iLC = 1:4%nLCs
%     for iImg = 1:nImg
%         tmpy = Lx(:,iLC);
%         tmpx = X(:,iImg);
%         r = corrcoef(tmpx,tmpy.');
%         LC_img_loadings(iImg,iLC) = r(1,2);
%         clear tmpy tmpx r
%     end
% end
% Dani: suggestion for faster computation of the identical matrix:
corr_Lx_X = corr(Lx,X)';
% Val: the 2 procedures give me different results (?!)

% Pearson's correlations between Ly and Y (behavior/design structure coefficients) 
for iLC = 1:nLCs
    for iB = 1:nBehav
        tmpy = Ly(:,iLC);
        tmpx = Y(:,iB);
        r = corrcoef(tmpx,tmpy.');
        LC_behav_loadings(iB,iLC) = r(1,2);
        clear tmpy tmpx r
    end
end
% Dani: suggestion for faster computation of the identical matrix:
% I believe that these correlations should be computed for each group
% separately as well, right?
% Val: yes, those using Ly

iter = 1;

for iG = 1:nGroups
    idx = iter:iter + nBehav - 1;
    this_groupID = find(grouping == groupIDs(iG));
    corr_Ly_Y(idx,:) = corr(Ly(this_groupID,:),Y(this_groupID,:))';
    iter = iter + nBehav;
end

% Pearson's correlations between Lx and Y
corr_Lx_Y = corr(Lx,Y);

% Pearson's correlations between Ly and X
%%%%% Adapt for groups because using Ly based on U which can be group-specific
corr_Ly_X = corr(Ly,X);
