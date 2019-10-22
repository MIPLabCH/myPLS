function Yp = myPLS_permuteY(Y,grouping,grouped_perm)

% This function permutes the Y matrix during permutation testing
%
% Inputs:
% - Y               : N x B vector, N is #subjects, B is #behavior/design
%                     scores, behavior/design data (normalized)
% - grouping        : N x 1 vector, subject grouping (e.g. diagnosis)
%                     e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                     subject 3 belongs to group 2
% - grouped_perm    : binary variable indicating if groups should be 
%                     considered during permutation of Y
%                     0 = permutations ignoring grouping
%                     1 = permutations within group
%
% Outputs:
% - Yp              : N x B vector, permuted behavior/design data

% Remove grouping information if it should not be considered for
% permutations
if ~grouped_perm
    grouping = ones(size(grouping));
end

% Get number of groups and their indices
groupIDs = unique(grouping);
nGroups = length(groupIDs);

Yp = zeros(size(Y));
for iG = 1:nGroups
    groupID = find(grouping == groupIDs(iG));   

    % Permute subjects of this group
    thisY = Y(groupID,:);
    perm_order = randperm(size(thisY,1));
    thisYp = thisY(perm_order,:);

    % Save permuted subjects in the right order
    Yp(groupID,:) = thisYp;
end
