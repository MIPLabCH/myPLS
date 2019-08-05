function X = myPLS_norm(X,grouping_PLS,mode)

% Normalize data
%
% IN:
%   X            : N x M matrix, N is #subjects, M is #imaging variables
%   grouping_PLS : N x 1 vector, subject grouping for PLS analysis
%                     e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                     subject 3 belongs to group 2.
%   mode         : normalization option 
%                   0 = no normalization
%                   1 = zscore across all subjects
%                   2 = zscore within groups (default)
%                   3 = std normalization across subjects (no centering)
%                   4 = std normalization within groups (no centering)
%
% OUT:
%   X : normalized X matrix
%
%
% function based on myPLS_norm by Dimitri Van De Ville
%   modifications D. Zoeller (Aug 2019): 
%       - adapted to work also for different group labels than 1,2,3,... 
%         (e.g. 2,3,4)
%       - removed group number input (can be derived from grouping vector)
%


% number and IDs of groups
groupIDs=unique(grouping_PLS);
nGroups=length(groupIDs);


switch mode
    case 1
        X=zscore(X);
    case 2
        for iG=1:nGroups
            idx=find(grouping_PLS==groupIDs(iG));
            X(idx,:)=zscore(X(idx,:));
        end
    case 3
        X2=sqrt(mean(X.^2,1));
        X=X./repmat(X2,[size(X,1) 1]);
    case 4
        for iG=1:nGroups
            idx=find(grouping_PLS==groupIDs(iG));
            X2=sqrt(mean(X(idx,:).^2,1));
            X(idx,:)=X(idx,:)./repmat(X2,[size(X(idx,:),1) 1]);
        end
end
