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
    sel_idx=iter:iter + nBehav - 1;
    Usel = U(sel_idx,:);        
    for iG2 = 1:nGroups
        groupID_find = find(grouping == groupIDs(iG2));
        Ysel = Y(groupID_find,:);
        Lyy(iG,groupID_find,:) = Ysel * Usel;
    end        
    iter = iter + nBehav;
end
for iG = 1:nGroups
    groupID_find = find(grouping == groupIDs(iG));
    Ly(groupID_find,:) = Lyy(iG,groupID_find,:);
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



% Contribution of original variables to LVs
% % Brain & behavior structure coefficients (Correlations imaging/behavior variables - brain/behavior scores)
% 
% clear myBrainStructCoeff myBehavStructCoeff
% 
% % Brain structure coefficients
% for iter_lv = 1:numSignifLVs
%     this_lv = mySignifLVs(iter_lv);
%     
%     for iter_img = 1:size(X,2)
%         clear tmpy tmpx r p
%         tmpx = X(:,iter_img);
%         tmpy = Lx(:,this_lv);
%         [r,p] = corrcoef(tmpx,tmpy.');
%         myBrainStructCoeff(iter_img,iter_lv) = r(1,2);
%     end
%     % Dani: same here, you could use the function corr() to compute
%     faster correlations
%     
% end
% 
% % Behavior structure coefficients
% for iter_lv = 1:numSignifLVs
%     this_lv = mySignifLVs(iter_lv);
% 
%     for iter_behav = 1:size(Y,2),
%         clear tmpy tmpx r p
%         tmpx = Y(:,iter_behav);
%         tmpy = Ly(:,this_lv);        
%         [r,p] = corrcoef(tmpx,tmpy.');
%         myBehavStructCoeff(iter_behav,iter_lv) = r(1,2);
%     end
%     % Dani: same here, you could use the function corr() to compute
%     faster correlations
% 
% end

