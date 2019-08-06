function Sp_vect = myPLS_permutations(X,Y,U,grouping,pls_opts)
%
% Permutation testing over singular values obtained with PLS
% Rows (subjects) of Y are permuted within each diagnostic group
%
% Inputs:
% - X          : N x M matrix, N is #subjects, M is #imaging
% - Y          : N x B matrix, B is #behaviors
% - U          : B x L matrix, L is #latent components, behavior saliences
% - grouping   : N x 1 vector, subject (diagnostic) grouping 
%                   e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                   subject 3 belongs to group 2.
% - pls_opts : options for the PLS analysis
%                necessary fields for this function:
%       - .nPerms : number of permutations
%       - .grouped_perm : binary variable indicating if groups should be 
%               considered during the permutations
%              0 = permutations ignoring grouping
%              1 = permutations within group
%       - .grouped_PLS         : binary variable indicating if groups
%                                should be considered when computing R
%              0 = PLS will computed over all subjects
%              1 = R will be constructed by concatenating group-wise
%                  covariance matrices ( as in conventional behavior PLS)
%
%
% Outputs:
% - Sp_vect             : matrix with bootstrapping salience values, from
%                         these, the p-values can be computed


% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

disp('... Permutations ...')
for iter_perm = 1:pls_opts.nPerms
    
    % Display number of permutations (every 50 samples) - (Dani: my personal
    % preference is to include a skipped line once in a while)
    if mod(iter_perm,20) == 0, fprintf('%d ',iter_perm); end
    if ~mod(iter_perm,200); fprintf('\n'); end
    
    % Leave X unchanged (no need to permute both X and Y matrices)
    Xp = X; % X is already normalized
    
    % Permute Y (within each diagnostic group, or across all subjects,
    % according to setup)
    Yp = myPLS_permuteY(Y,grouping,pls_opts.grouped_perm);
    
    % Cross-covariance matrix between X and permuted Y
    Rp = myPLS_cov(Xp,Yp,grouping,pls_opts.grouped_PLS);
    
    % Singular value decomposition of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat;
    Sp = sqrt(sum(Up.^2));
    
    % Keep singular values for sample distribution of singular values
    Sp_vect(:,iter_perm) = Sp';    
end

disp(' ')
