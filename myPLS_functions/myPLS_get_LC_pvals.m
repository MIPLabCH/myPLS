function LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts)

% This function computed pvalues for all latent components (LCs) using the
% permuted singular values 
%
% Inputs
% - Sp_vect    : L x #permutations matrix, L is #LCs, permuted singular values, 
%                used to compute p-values to assess LCs' significance                         
% - S          : L x L matrix, singular values (diagonal matrix)
% - pls_opts   : options for the PLS analysis
%       - .nPerms : number of permutations
%
% Outputs:
% - LC_pvals   : L x 1 vector, p-values for all LCs

% Number of permutations with S_perm greater than S
S_mat = repmat(diag(S),1,pls_opts.nPerms);
sp = sum(Sp_vect >= S_mat,2);

% Compute p-values for each LC - approximation by counting the number of
% permuted singular values above the measured singular value
LC_pvals = (sp + 1) ./ (pls_opts.nPerms + 1);

signif_LC = find(LC_pvals<0.05); % index of significant LVs
nSignifLC = size(signif_LC,1); % number of significant LVs

% Display significant LCs
disp([num2str(nSignifLC) ' significant LC(s) - p < 0.05 (uncorrected!)']);
for iLC = 1:nSignifLC
    this_lc = signif_LC(iLC);
    disp(['LC' num2str(this_lc) ' - p=' num2str(LC_pvals(this_lc),'%0.3f') ]);
end

disp(' ')
