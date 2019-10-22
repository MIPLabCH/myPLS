function Sp_vect = myPLS_permutations(X,Y,U,grouping,pls_opts)

% This function computed permutation testing over singular values obtained
% with PLS. Rows (subjects) of Y are permuted in each iteration.
%
% Inputs:
% - X          : N x M matrix, N is #subjects, M is #imaging variables, 
%                imaging data (normalized)
% - Y          : N x B matrix, B is #behavior scores, behavior/design data
%                (normalized)
% - U          : B x L matrix, L is #latent components, behavior/design saliences
% - grouping   : N x 1 vector, subject grouping (e.g. diagnosis)
%                e.g. [1,1,2] = subjects 1 and 2 belong to group 1,
%                subject 3 belongs to group 2
% - pls_opts   : options for the PLS analysis
%                Necessary fields for this function:
%       - .nPerms       : number of permutations
%       - .grouped_perm : binary variable indicating if groups should be 
%                considered during the permutations
%                0 = permutations ignoring grouping
%                1 = permutations within group
%       - .grouped_PLS  : binary variable indicating if groups
%                should be considered when computing R
%                0 = PLS will computed over all subjects
%                1 = R will be constructed by concatenating group-wise
%                covariance matrices ( as in conventional behavior PLS)
%
% Outputs:
% - Sp_vect    : L x #permutations matrix, permuted singular values, used to compute
%                p-values to assess LCs' significance                         

% Set up random number generator
rng(1);

% Check that dimensions of X & Y are correct
if(size(X,1) ~= size(Y,1))
    error('Input arguments X and Y should have the same number of rows');
end

disp('... Permutations ...')
for iP = 1:pls_opts.nPerms
        
    % Leave X unchanged (no need to permute both X and Y matrices)
    Xp = X; % X is already normalized
    
    % Permute Y (within each group or across all subjects)
    Yp = myPLS_permuteY(Y,grouping,pls_opts.grouped_perm);
    
    % Generate cross-covariance matrix between X and permuted Y
    Rp = myPLS_cov(Xp,Yp,grouping,pls_opts.grouped_PLS);
    
    % Singular value decomposition of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    
    % Procrustas transform (correction for axis rotation/reflection)
    rotatemat = rri_bootprocrust(U, Up);
    Up = Up * Sp * rotatemat;
    Sp = sqrt(sum(Up.^2));
    
    % Keep singular values for sample distribution of singular values
    Sp_vect(:,iP) = Sp'; 
    
end

if mod(iP,200); fprintf('\n'); end
disp(' ')
