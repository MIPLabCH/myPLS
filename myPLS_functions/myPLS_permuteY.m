function Yp = myPLS_permuteY(Y,grouping,consider_groups)

% remove grouping information if it should not be considered for
% permutations
if ~consider_groups
    grouping=ones(size(grouping));
end

% get number of groups and their indices
groupIDs=unique(grouping);
nGroups=length(groupIDs);

Yp = zeros(size(Y));
for iG = 1:nGroups
    groupID=grouping == groupIDs(iG);

    % permute subjects of this group
    thisY = Y(groupID,:);
    perm_order = randperm(size(thisY,1));
    thisYp = thisY(perm_order,:);

    % save permuted subjects in the right places
    Yp(groupID,:) = thisYp;
end