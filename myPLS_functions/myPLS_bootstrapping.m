function boot_results = myPLS_bootstrapping(X0,Y0,U,V,S,grouping,pls_opts)

% This function computes bootstrap resampling with replacement on X and Y, 
% with option to either ignore or account for groups (e.g. diagnosis).
% The PLS analysis is re-computed for each bootstrap sample. 
%
% Inputs:
% - X0             : N x M matrix, N is #subjects, M is #imaging,
%                    brain data (not normalized!)
% - Y0             : N x B matrix, B is #behaviors/design scores,
%                    behavior/design data (not normalized!)
% - U              : B x L matrix, L is #latent components (LCs), 
%                    behavior/design saliences
% - V              : B x L matrix, imaging saliences
% - S              : L x L matrix, singular values (diagonal matrix)
% - grouping       : N x 1 vector, subject grouping (e.g. diagnosis)
%                    e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                    subject 3 belongs to group 2
% - pls_opts       : options for the PLS analysis
%                    Necessary fields for this function:
%       - .nBootstraps : number of bootstrap samples
%       - .grouped_PLS : binary variable indicating if groups
%                        should be considered when computing R
%              0 = PLS will computed over all subjects
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices (as in conventional behavior PLS)
%       - .grouped_boot : binary variable indicating if groups should be 
%                         considered during bootstrapping
%              0 = bootstrapping ignoring grouping
%              1 = bootstrapping within group
%       - .boot_procrustes_mod : mode for bootstrapping procrustes transform
%              1 = standard (rotation computed only on U)
%              2 = average rotation of U and V
%       - .save_boot_resampling: indicator whether to save bootstrap
%                                resampling data or not
%              0 = no saving of bootstrapping resampling data
%              1 = save bootstrapping resampling data
%       - .normalization_img    : normalization options for imaging data
%       - .normalization_behav  : normalization options for behavior/design data
%              0 = no normalization
%              1 = zscore across all subjects
%              2 = zscore within groups (default)
%              3 = std normalization (no centering) across subjects 
%              4 = std normalization (no centering) within groups
%
% Outputs:
% - boot_results : struct containing all results from bootstrapping
%       - .Ub_vect    : B x S x P matrix, S is #LCs, P is #bootstrap samples,
%                       bootstrapped behavior saliences for all LCs
%       - .Vb_vect    : M x S xP matrix, bootstrapped brain saliences for all LCs
%       - .Lxb,.Lyb,.LC_img_loadings_boot,.LC_behav_loadings_boot :
%                       bootstrapping scores (see myPLS_get_PLSscores for
%                       details) 
%       - .*_mean : mean of bootstrapping distributions
%       - .*_std : standard deviation of bootstrapping distributions
%       - .*_lB : lower bound of 95% confidence interval of bootstrapping distributions
%       - .*_uB : upper bound of 95% confidence interval of bootstrapping distributions

% Set up random number generator
rng(1);

% Check that dimensions of X & Y are correct
if(size(X0,1) ~= size(Y0,1))
    error('Input arguments X and Y should have the same number of rows');
end

nBehav = size(Y0,2);
nGroups = size(unique(grouping),1);


disp('... Bootstrapping ...')

% Compute the bootstrapping orders (under consideration of the grouping, if asked for)
all_boot_orders = myPLS_get_boot_orders(pls_opts.nBootstraps,grouping,pls_opts.grouped_boot);

% Run PLS in each bootstrap sample
for iB = 1:pls_opts.nBootstraps
    
    % Resampling X
    Xb = X0(all_boot_orders(:,iB),:);
    Xb = myPLS_norm(Xb,grouping,pls_opts.normalization_img);
    
    % Resampling Y
    Yb = Y0(all_boot_orders(:,iB),:);
    Yb = myPLS_norm(Yb,grouping,pls_opts.normalization_behav);
    
    % Generate cross-covariance matrix between resampled X and Y
    Rb = myPLS_cov(Xb,Yb,grouping,pls_opts.grouped_PLS);
    
    % Singular value decomposition of Rb
    [Ub,Sb,Vb] = svd(Rb,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    switch pls_opts.boot_procrustes_mod
        case 1
            % Computed on U only
            rotatemat_full = rri_bootprocrust(U, Ub);
            
            % Rotate and re-scale Ub and Vb
            Vb = Vb * Sb * rotatemat_full;
            Ub = Ub * Sb * rotatemat_full;
            
            Vb = Vb./repmat(diag(S)',size(Vb,1),1);
            Ub = Ub./repmat(diag(S)',size(Ub,1),1);
        case 2
            % Computed on both U and V
            rotatemat1 = rri_bootprocrust(U, Ub);
            rotatemat2 = rri_bootprocrust(V, Vb);
    
            % Full rotation
            rotatemat_full = (rotatemat1 + rotatemat2)/2;
            
            % Apply full rotation to Vb and Ub
            Vb = Vb * rotatemat_full;
            Ub = Ub * rotatemat_full;
            
        otherwise
            error('invalid value in pls_opts.boot_procrustes_mod!');
    end
    
    % Vectors with all bootstrap samples -> needed for percentile computation
    Ub_vect(:,:,iB) = Ub;
    Vb_vect(:,:,iB) = Vb;
    
    % Compute bootstrapping PLS scores
    [Lxb,Lyb,corr_Lxb_Xb,corr_Lyb_Yb,corr_Lxb_Yb,corr_Lyb_Xb] = ...
        myPLS_get_PLS_scores_loadings(Xb,Yb,Vb,Ub,grouping,pls_opts);
    
    boot_results.Lxb(:,:,iB) = Lxb;
    boot_results.Lyb(:,:,iB) = Lyb;
    boot_results.LC_img_loadings_boot(:,:,iB) = corr_Lxb_Xb;
    boot_results.LC_behav_loadings_boot(:,:,iB) = corr_Lyb_Yb;
    
end

% Compute bootstrapping statistics
boot_stats = myPLS_bootstrap_stats(Ub_vect,Vb_vect,boot_results);

% Save all the statistics fields in the boot_results (I am coding it like
% this to facilitate adding more stats
fN = fieldnames(boot_stats);
for iF = 1:length(fN)
    boot_results.(fN{iF}) = boot_stats.(fN{iF});
end

% Save bootstrapping resampling data if asked for
if pls_opts.save_boot_resampling
    boot_results.Ub_vect = Ub_vect;
    boot_results.Vb_vect = Vb_vect;
else
    boot_results = rmfield(boot_results,'LC_img_loadings_boot');
    boot_results = rmfield(boot_results,'LC_behav_loadings_boot');
end


if mod(iB,200); fprintf('\n'); end
disp(' ')
