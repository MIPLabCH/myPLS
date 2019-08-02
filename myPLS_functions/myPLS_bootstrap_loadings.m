function [LC_img_loadings_boot,LC_behav_loadings_boot] = myPLS_bootstrap_loadings(X,Y,U,V,grouping,signif_LC,nBootstraps,nGroups_PLS,grouping_PLS,normalization_img,normalization_behav)
% 
% This function computes bootstrap resampling with replacement on X and Y, 
% accounting for diagnostic groups.
% The PLS analysis is re-computed for each bootstrap sample. 
% The imaging & behavior loadings are also computed for each bootstrap sample
% for all significant latent components (LCs). 
%
% Inputs:
% - X                    : N x M matrix, N is #subjects, M is #imaging
% - Y                    : N x B matrix, B is #behaviors
% - U                    : B x L matrix, L is #latent components, behavior saliences
% - V                    : B x L matrix, imaging saliences
% - grouping             : N x 1 vector, subject (diagnostic) grouping 
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
% - signif_LC            : significant latent components to plot
% - nBootstraps          : number of bootstrap samples
% - nGroups_PLS          : number of groups for PLS analysis, will determine
%                          if each group has its own set of saliences
%                          1 = the cross-covariance matrix is computed across all subjects
%                          2 = the cross-covariance matrix is computed within each group, and
%                          each group has its set of behavior saliences 
% - grouping_PLS         : N x 1 vector, subject grouping for PLS analysis
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
% - normalization_img    : normalization options for FC data
% - normalization_behav  : normalization options for behavior data
%                          0 = no normalization
%                          1 = zscore across all subjects
%                          2 = zscore within groups (default)
%                          3 = std normalization across subjects (no centering)
%                          4 = std normalization within groups (no centering)
%
% Outputs:
% - LC_img_loadings_boot   : M x S x P matrix, S is #significant LCs, P is
%                            #bootstrap samples, bootstrapped imaging loadings
%                            for significant LCs
% - LC_behav_loadings_boot : B x S x P matrix, bootstrapped behavior loadings
%                            for significant LCs

% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

nBehav = size(Y,2);
nGroups = size(unique(grouping),1);

% Bootstrap with replacement within each diagnostic group
all_boot_orders = [];
for iter_group = 1:nGroups
    nSubj_group = size(find(grouping==iter_group),1);
    [boot_order,~] = rri_boot_order(nSubj_group,1,nBootstraps);
    all_boot_orders = [all_boot_orders; boot_order];
    clear boot_order nSubj_group
end


for iter_boot = 1:nBootstraps
    
    if mod(iter_boot,50)==0, disp(num2str(iter_boot)); end
    
    % Resampling X
    Xb = X(all_boot_orders(:,iter_boot),:);
    Xb = myPLS_norm(Xb,nGroups_PLS,grouping_PLS,normalization_img);
    
    % Resampling Y
    Yb = Y(all_boot_orders(:,iter_boot),:);
    Yb = myPLS_norm(Yb,nGroups_PLS,grouping_PLS,normalization_behav);
    
    % Cross-covariance matrix between resampled X and Y
    Rb = myPLS_cov(Xb,Yb,nGroups_PLS,grouping_PLS);
    
    % Singular value decomposition of Rp
    [Ub,~,Vb] = svd(Rb,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    % Computed on both U and V to
    rotatemat = rri_bootprocrust(U, Ub);
    rotatemat2 = rri_bootprocrust(V, Vb);
    
    % Full rotation
    rotatemat_full = (rotatemat + rotatemat2)/2;
    Vb = Vb * rotatemat_full;
    Ub = Ub * rotatemat_full;
    
    % Imaging & behavior composite scores
    Lxb = Xb * Vb;
    
    if nGroups_PLS == 1
        
        Lyb = Yb * Ub;
        
    elseif nGroups_PLS == 2
        Lyb = nan(size(Lxb));
        
        iter = 1;
        for iter_group = 1:nGroups_PLS
            Ubsel = Ub(iter:iter + nBehav - 1,:);
            for iter_group2 = 1:nGroups_PLS
                idx = find(grouping_PLS == iter_group2);
                Ybsel = Yb(idx,:);
                Lyyb(iter_group,idx,:) = Ybsel * Ubsel;
            end
            iter = iter + nBehav;
        end
        
        for iter_group = 1:nGroups_PLS
            idx = find(grouping_PLS == iter_group);
            first = idx(1);
            last = idx(end);
            Lyb(first:last,:) = Lyyb(iter_group,first:last,:);
        end
        
        clear Lyy Usel Ysel idx iter first last
    end
    
    % Loadings
    for iter_lc = 1:size(signif_LC,1)
        this_lc = signif_LC(iter_lc);
        
        % Imaging
        for iter_img = 1:size(Xb,2)
            tmpy = Lxb(:,this_lc);
            tmpx = Xb(:,iter_img);
            [r,~] = corrcoef(tmpx,tmpy.');
            these_img_loadings(iter_img,iter_lc) = r(1,2);
            clear tmpy tmpx r
        end
        
        % Behavior
        for iter_behav = 1:size(Y,2)
            tmpy = Lyb(:,this_lc);
            tmpx = Yb(:,iter_behav);
            [r,~] = corrcoef(tmpx,tmpy.');
            these_behav_loadings(iter_behav,iter_lc) = r(1,2);
            clear tmpy tmpx r
        end
    end
    
    LC_img_loadings_boot(:,:,iter_boot) = these_img_loadings;
    LC_behav_loadings_boot(:,:,iter_boot) = these_behav_loadings;
    
end
