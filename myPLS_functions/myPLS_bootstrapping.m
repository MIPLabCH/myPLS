function [Ub_vect,Vb_vect] = myPLS_bootstrapping(X0,Y0,U,V,grouping,pls_opts)

% 
% This function computes bootstrap resampling with replacement on X and Y, 
% with option to either ignore or accounte for diagnostic groups.
% The PLS analysis is re-computed for each bootstrap sample. 
%
% Inputs:
% - X0                   : N x M matrix, N is #subjects, M is #imaging, brain data not normalized!
% - Y0                   : N x B matrix, B is #behaviors/design-scores, behavior/design data not normalized!
% - U                    : B x L matrix, L is #latent components, behavior/design saliences
% - V                    : B x L matrix, imaging saliences
% - grouping             : N x 1 vector, subject (diagnostic) grouping 
%                          e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                          subject 3 belongs to group 2.
% - pls_opts : options for the PLS analysis
%                necessary fields for this function:
%       - .nBootstraps          : number of bootstrap samples
%       - .grouped_PLS         : binary variable indicating if groups
%                                should be considered when computing R
%              0 = PLS will computed over all subjects
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices ( as in conventional behavior PLS)
%       - .grouped_boot : binary variable indicating if groups should be 
%               considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - .boot_procrustes_mod : mode for bootstrapping procrustes transform
%              1 = standard (rotation computed only on U)
%              2 = average rotation of U and V
%       - .normalization_img    : normalization options for FC data
%       - .normalization_behav  : normalization options for behavior data
%              0 = no normalization
%              1 = zscore across all subjects
%              2 = zscore within groups (default)
%              3 = std normalization across subjects (no centering)
%              4 = std normalization within groups (no centering)
%
% Outputs:
% - Ub_vect : P x B x S matrix, S is #LCs, P is #bootstrap samples,
%             bootstrapped behavior saliences for all LCs
% - Vb_vect : P x M x S matrix, bootstrapped brain saliences for all LCs



% Check that dimensions of X & Y are correct
if(size(X0,1) ~= size(Y0,1))
    error('Input arguments X and Y should have the same number of rows');
end

nBehav = size(Y0,2);
nGroups = size(unique(grouping),1);


disp('... Bootstrapping ...')

% compute the bootstrapping orders (under consideration of the grouping, if asked for)
all_boot_orders = myPLS_get_boot_orders(pls_opts.nBootstraps,grouping,pls_opts.grouped_boot);


for iter_boot = 1:pls_opts.nBootstraps
    
    % Display number of Bootstraps (every 50 samples) - (Dani: my personal
    % preference is to include a skipped line once in a while)
    if mod(iter_boot,20) == 0, fprintf('%d ',iter_boot); end
    if ~mod(iter_boot,200); fprintf('\n'); end
    
    % Resampling X
    Xb = X0(all_boot_orders(:,iter_boot),:);
    Xb = myPLS_norm(Xb,grouping,pls_opts.normalization_img);
    
    % Resampling Y
    Yb = Y0(all_boot_orders(:,iter_boot),:);
    Yb = myPLS_norm(Yb,grouping,pls_opts.normalization_behav);
    
    % Cross-covariance matrix between resampled X and Y
    Rb = myPLS_cov(Xb,Yb,grouping,pls_opts.grouped_PLS);
    
    % Singular value decomposition of Rp
    [Ub,~,Vb] = svd(Rb,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    switch pls_opts.boot_procrustes_mod
        case 1
            % Computed on only U
            rotatemat_full = rri_bootprocrust(U, Ub);
        case 2
            % Computed on both U and V
            rotatemat1 = rri_bootprocrust(U, Ub);
            rotatemat2 = rri_bootprocrust(V, Vb);
    
            % Full rotation
            rotatemat_full = (rotatemat1 + rotatemat2)/2;
        otherwise
            error('invalid value in pls_opts.boot_procrustes_mod!');
    end
    
    Vb = Vb * rotatemat_full;
    Ub = Ub * rotatemat_full;
    
    % vectors with all bootstrap samples --> needed for percentile
    % computation
    Ub_vect(iter_boot,:,:) = Ub;
    Vb_vect(iter_boot,:,:) = Vb;
end


disp(' ')
