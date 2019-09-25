function all_boot_orders = myPLS_get_boot_orders(nBootstraps,grouping,grouped_boot)

% This function computes bootstrapping sampling orders for bootstrapping 
% with resampling. Note that if groups are considered for resampling, 
% groups can be not ordered.
%
% Inputs:
% - nBootstraps     : #bootstrap samples
% - grouping        : N x 1 vector, subject grouping (e.g., diagnosis)
% - grouped_boot    : binary variable indicating if groups should be 
%                     considered during bootstrapping
%                     0 = bootstrapping ignoring grouping
%                     1 = bootstrapping within group
%
% Outputs:
% - all_boot_orders : N x bootstrap samples, subject index in bootstrap
%                     samples


% Ignore grouping information if asked for
if ~grouped_boot
    grouping = ones(size(grouping));
end

% Get number of groups and their indices
groupIDs = unique(grouping);
nGroups = length(groupIDs);

% Number of subjects
nSubj = length(grouping); 

% Get bootstrap subject sampling
all_boot_orders = nan(nSubj,nBootstraps);

for iG = 1:nGroups
    % Get the IDs of the subjects in this group
    groupID = find(grouping == groupIDs(iG));
    num_subj_group = length(groupID);

    % Compute bootstrapping orders for the number of subjects in this group
    [boot_order_tmp,~] = rri_boot_order(num_subj_group,1,nBootstraps);

    % Change the order of the group indices according to the bootstrapping
    all_boot_orders(groupID,:) = groupID(boot_order_tmp);
end
