function R = myPLS_cov(X,Y,grouping,grouped_PLS)

% This function generates a cross-covariance matrix (stacked) for the PLS
% analysis
% 
% Inputs:
% - X            : N x M matrix, N is #subjects, M is #imaging variables
% - Y            : N x B matrix, B is #behaviors
% - grouping     : N x 1 vector, subject grouping (e.g. diagnosis)
%                  e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                  subject 3 belongs to group 2
% - grouped_PLS  : binary variable indicating if groups should be 
%                  considered when computing R
%              0 = PLS will computed over all subjects (ignoring grouping)
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices (as in conventional behavior PLS)
%
% Outputs:
% - R            : cross-covariance matrix
%
%
% Function based on myPLS_cov by Dimitri Van De Ville
% Modifications by D. Zoeller (Aug 2019): 
%   - adapted to work also for different group labels than 1,2,3,... (e.g. 2,3,4)
%   - removed group number input (can be derived from grouping vector)


if ~grouped_PLS
    grouping = ones(size(grouping));
end

% Number and IDs of groups
groupIDs = unique(grouping);
nGroups = length(groupIDs);

for iG = 1:nGroups
    Ysel = Y(grouping == groupIDs(iG),:);
    Xsel = X(grouping == groupIDs(iG),:);
    
    R0 = Ysel.'*Xsel;

    if ~exist('R')
        R = R0;
    else
        R = [R; R0];
    end
end
