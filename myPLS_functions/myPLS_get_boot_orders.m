function all_boot_orders = myPLS_get_boot_orders(nBootstraps,grouping,grouped_boot)

% function to compute bootstrapping sampling orders for bootstrapping with
% resampling
% if groups are considered for resampling, groups can be not ordered


% ignoring grouping information if asked for
if ~grouped_boot
    grouping=ones(size(grouping));
end


% get number of groups and their indices
groupIDs=unique(grouping);
nGroups=length(groupIDs);

% number of subjects
nSubj = length(grouping); 

% Get bootstrap subject sampling
all_boot_orders = nan(nSubj,nBootstraps);
for iG = 1:nGroups
    % getting the IDs of the subjects in this group
    groupID_find=find(grouping == groupIDs(iG));
    num_subj_group = length(groupID_find);

    % compute bootstrapping orders for the number of subjects in this group
    [boot_order_tmp,~] = rri_boot_order(num_subj_group,1,nBootstraps);

    % change the order of the group indices according to the bootstrapping
    all_boot_orders(groupID_find,:) = groupID_find(boot_order_tmp);
end

