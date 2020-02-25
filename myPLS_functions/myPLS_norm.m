function [X, meanX, stdX] = myPLS_norm(X,grouping,mode)

% This function applies normalization on data before running PLS
%
% Inputs:
% - X            : N x V matrix, N is #subjects, V is #variables
% - grouping     : N x 1 vector, subject grouping (e.g. diagnosis)
%                  e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                  subject 3 belongs to group 2
% - mode         : normalization option 
%                  0 = no normalization
%                  1 = zscore across all subjects
%                  2 = zscore within groups (default)
%                  3 = std normalization (no centering) across subjects 
%                  4 = std normalization (no centering) within groups 
%
% Outputs:
% - X            : normalized X matrix
%
%
% Function based on myPLS_norm by D. Van De Ville
% Modifications by D. Zoeller (Aug 2019): 
%   - adapted to work also for different group labels than 1,2,3,... (e.g. 2,3,4)
%   - removed group number input (can be derived from grouping vector)


% Number and IDs of groups
groupIDs = unique(grouping);
nGroups = length(groupIDs);


switch mode
    case 1
        meanX=mean(X);
        stdX=std(X);
        X = zscore(X);
    
    case 2
        for iG = 1:nGroups
            idx = find(grouping == groupIDs(iG));
            meanX(iG,:)=mean(X(idx,:));
            stdX(iG,:)=std(X(idx,:));
            X(idx,:) = zscore(X(idx,:));
        end
        
    case 3
        meanX=mean(X);
        stdX=sqrt(mean(X.^2,1));
        X2 = stdX;
        X = X./repmat(X2,[size(X,1) 1]);
    
    case 4
        for iG=1:nGroups
            idx = find(grouping == groupIDs(iG));
            meanX(iG,:)=mean(X(idx,:));
            stdX(iG,:)=sqrt(mean(X(idx,:).^2,1));
            X2 = stdX(iG,:);
            X(idx,:) = X(idx,:)./repmat(X2,[size(X(idx,:),1) 1]);
        end
end
