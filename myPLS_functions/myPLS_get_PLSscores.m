function [Lx,Ly,LC_img_loadings,LC_behav_loadings] = myPLS_get_PLSscores(X,Y,V,U,grouping,pls_opts)

% Function to compute PLS scores
% I propose to simply compute all scores that we can imagine, and then only
% give options which ones to plot in the end


if ~pls_opts.grouped_PLS
    grouping=ones(size(grouping));
end

% number and IDs of groups
groupIDs=unique(grouping);
nGroups=length(groupIDs);

% number of behavior scores
nBehav=size(Y,2);

% Number of latent components (LC)
nLCs = size(V,2); 

% number of imaging measures
nImg = size(V,1);  


%% copied (and slightly modified) from original myPLS_analysis function:

% Imaging & behavioral composite scores
Lx = X * V;

% Dani: I removed the selection depending on the number of groups, I believe
% it should work like this even for one group
iter = 1;
for iG = 1:nGroups        
    Usel = U(iter:iter + nBehav - 1,:);        
    for iG2 = 1:nGroups
        groupID_find = find(grouping == iG2);
        Ysel = Y(groupID_find,:);
        Ly(groupID_find,:) = Ysel * Usel;
    end        
    iter = iter + nBehav;
end


% Imaging loadings (Pearson's correlations between Lx and X)
% for iLC = 1:nLCs
%     for iter_img = 1:nImg
%         tmpy = Lx(:,iLC);
%         tmpx = X(:,iter_img);
%         r = corrcoef(tmpx,tmpy.');
%         LC_img_loadings(iter_img,iLC) = r(1,2);
%         clear tmpy tmpx r
%     end
% end
% Dani: suggestion for faster computation of the identical matrix:
LC_img_loadings=corr(Lx,X)';

% Behavior loadings (Pearson's correlations between Ly and Y)
% for iLC = 1:nLCs
%     for iter_behav = 1:nBehav
%         tmpy = Ly(:,iLC);
%         tmpx = Y(:,iter_behav);
%         r = corrcoef(tmpx,tmpy.');
%         LC_behav_loadings(iter_behav,iLC) = r(1,2);
%         clear tmpy tmpx r
%     end
% end
% Dani: suggestion for faster computation of the identical matrix:
LC_behav_loadings=corr(Ly,Y)';


