function boot_stats = myPLS_bootstrap_stats(Ub_vect,Vb_vect,boot_results)

% This function computes the mean, standard deviation and confidence intervals
% across the imaging & behavior/design loadings obtained with bootstrap 
% resampling. Bootstrap scores are not computed at this stage, as they 
% depend on the PLS results.
%
% Inputs:
% - Ub_vect                 : B x L x P matrix, B is #behavior/design scores, L is #LCs,
%                             P is #bootstrap samples, bootstrapped behavior/design saliences
% - Vb_vect                 : M x S x P matrix, bootstrapped imaging saliences
% - boot_results            : struct containing results from bootstrapping
%     - .Lxb,.Lyb,.LC_img_loadings_boot,.LC_behav_loadings_boot :
%                bootstrapping scores (see myPLS_get_PLS_scores_loadings for details)
% - LC_behav_loadings_boot  : B x S x P matrix, S is #LCs, P is # bootstrap samples,
%                             bootstrapped behavior/design loadings
% - LC_img_loadings_boot    : M x S x P matrix, bootstrapped imaging loadings
%
% Outputs:
% - boot_stats : struct containing all computed bootstrapping measures
%     - .*_mean : mean of bootstrapping distributions
%     - .*_std : standard deviation of bootstrapping distributions
%     - .*_lB : lower bound of 95% confidence interval of bootstrapping distributions
%     - .*_uB : upper bound of 95% confidence interval of bootstrapping distributions

nBehav = size(boot_results.LC_behav_loadings_boot,1);
nImaging = size(boot_results.LC_img_loadings_boot,1); 
nBootstraps = size(boot_results.LC_behav_loadings_boot,3);

%% Behavior/design and imaging saliences

% Computing mean, std and CIs for bootstrap saliences
boot_stats.Ub_mean = mean(Ub_vect,3); % mean
boot_stats.Ub_std = std(Ub_vect,[],3); % standard deviation
boot_stats.Ub_lB = prctile(Ub_vect,2.5,3); % lower bound of 95% CI
boot_stats.Ub_uB = prctile(Ub_vect,97.5,3); % upper bound of 95% CI

boot_stats.Vb_mean = mean(Vb_vect,3); % mean
boot_stats.Vb_std = std(Vb_vect,[],3); % standard deviation
boot_stats.Vb_lB = prctile(Vb_vect,2.5,3); % lower bound of 95% CI
boot_stats.Vb_uB = prctile(Vb_vect,97.5,3); % upper bound of 95% CI

%% Behavior and imaging loadings

% Computing mean, std and CIs for bootstrap loadings
boot_stats.LC_behav_loadings_mean = mean(boot_results.LC_behav_loadings_boot,3); % mean
boot_stats.LC_behav_loadings_std = std(boot_results.LC_behav_loadings_boot,[],3); % standard deviation
boot_stats.LC_behav_loadings_lB = prctile(boot_results.LC_behav_loadings_boot,2.5,3); % lower bound of 95% CI
boot_stats.LC_behav_loadings_uB = prctile(boot_results.LC_behav_loadings_boot,97.5,3); % upper bound of 95% CI

boot_stats.LC_img_loadings_mean = mean(boot_results.LC_img_loadings_boot,3); % mean
boot_stats.LC_img_loadings_std = std(boot_results.LC_img_loadings_boot,[],3); % standard deviation
boot_stats.LC_img_loadings_lB = prctile(boot_results.LC_img_loadings_boot,2.5,3); % lower bound of 95% CI
boot_stats.LC_img_loadings_uB = prctile(boot_results.LC_img_loadings_boot,97.5,3); % upper bound of 95% CI
