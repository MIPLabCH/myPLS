function R = myPLS_cov(X,Y,subj_grouping)

% Compute PLS cross-covariance matrix (stacked)

% IN:
%   X            : N x M matrix, N is #subjects, M is #imaging variables
%   Y            : N x B matrix, B is #behaviors
%   grouping_PLS : N x 1 vector, subject grouping for PLS analysis
%                     e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                     subject 3 belongs to group 2.
%
% OUT:
%   R : cross-covariance matrix
%
%
% function based on myPLS_norm by Dimitri Van De Ville
%   modifications D. Zoeller (Aug 2019): 
%       - adapted to work also for different group labels than 1,2,3,... 
%         (e.g. 2,3,4)
%       - removed group number input (can be derived from grouping vector)



% number and IDs of groups
groupIDs=unique(subj_grouping);
nGroups=length(groupIDs);

for iG=1:nGroups
    Ysel = Y(subj_grouping==groupIDs(iG),:);
    Xsel = X(subj_grouping==groupIDs(iG),:);
    
    R0 = Ysel.'*Xsel;

    if ~exist('R')
        R = R0;
    else
        R = [R; R0];
    end
end

