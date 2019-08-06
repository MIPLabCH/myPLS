function LC_pvals = myPLS_get_LC_pvals(Sp_vect,S,pls_opts)

% number of permutations with S_perm greater than S
S_mat=repmat(diag(S),1,pls_opts.nPerms);
sp = sum(Sp_vect >= S_mat,2);

% compute the p-value of each LV - approximation by counting the number of
% permuted singular values above the measured singular value
LC_pvals = (sp + 1) ./ (pls_opts.nPerms + 1);

mySignifLVs = find(LC_pvals<0.05); % index of significant LVs
numSignifLVs = size(mySignifLVs,1); % number of significant LVs

% Display significant LVs
disp([num2str(numSignifLVs) ' significant LV(s) - p<0.05 (uncorrected!)']);
for iter_lv = 1:numSignifLVs
    this_lv = mySignifLVs(iter_lv);
    disp(['LV' num2str(this_lv) ' - p=' num2str(LC_pvals(this_lv),'%0.3f') ]);
end


disp(' ')
