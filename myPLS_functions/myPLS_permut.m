function pvals_LC = myPLS_permut(X,Y,U,S,nPerms,grouping,normalization_img,normalization_behav,nGroups_PLS,grouping_PLS)
%
% Permutation testing over singular values obtained with PLS
% Rows (subjects) of Y are permuted within each diagnostic group
%
% Inputs:
% - X                    : N x M matrix, N is #subjects, M is #imaging
% - Y                    : N x B matrix, B is #behaviors
% - U                    : B x L matrix, L is #latent components, behavior saliences
% - S                    : L x L matrix, singular values
% - nPerms               : number of permutations
% - grouping             : N x 1 vector, subject (diagnostic) grouping 
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
% - normalization_img    : normalization options for imaging data
% - normalization_behav  : normalization options for behavior data
%                          0 = no normalization
%                          1 = zscore across all subjects
%                          2 = zscore within groups (default)
%                          3 = std normalization across subjects (no centering)
%                          4 = std normalization within groups (no centering)
% - nGroups_PLS          : number of groups for PLS analysis, will determine
%                          if each group has its own set of saliences
%                          1 = the cross-covariance matrix is computed across all subjects
%                          2 = the cross-covariance matrix is computed within each group, and
%                          each group has its set of behavior saliences 
% - grouping_PLS         : N x 1 vector, subject grouping for PLS analysis
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
%
% Outputs:
% - pvals_LC             : p-values for each latent component


% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nSubj = size(X,1);
nGroups = size(unique(grouping),1);

grouping_PLS = ones(nSubj,1);

disp('... Permutations ...')
for iter_perm = 1:nPerms
    
    % Display number of permutations (every 50 permuts)
    if mod(iter_perm,50)==0, disp(num2str(iter_perm)); end
    
    % Normalization of X
    Xp = myPLS_norm(X,nGroups_PLS,grouping_PLS,normalization_img);
    
    % Permute Y within each diagnostic group
    Yp = [];
    for iter_group = 1:nGroups
        clear thisY thisYp
        thisY = Y(find(grouping==iter_group),:);
        perm_order = randperm(size(thisY,1));
        thisYp = thisY(perm_order,:);
        Yp = [Yp; thisYp];
    end
    
    % Normalization of Y
    Yp = myPLS_norm(Yp,nGroups_PLS,grouping_PLS,normalization_img);
    
    % Cross-covariance matrix between X and permuted Y
    Rp = myPLS_cov(Xp,Yp,nGroups_PLS,grouping_PLS);
    
    % Singular value decomposition of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat;
    Sp = sqrt(sum(Up.^2));
    
    % Keep singular values for sample distribution of singular values
    permsamp(:,iter_perm) = Sp';
    
    if iter_perm == 1
        sp= (Sp' >= diag(S));
    else
        sp = sp + (Sp'>= diag(S));
    end
    
end

% Compute p-values for each LC
pvals_LC = (sp + 1) ./ (nPerms + 1);