function boot_stats = myPLS_bootstrap_stats(Ub_vect,Vb_vect,boot_results)

% This function computes the mean, std deviation and confidence intervals
% across the imaging & behavior/design loadings obtained with bootstrap 
% resampling (myPLS_bootstrapping).
% Bootstrap scores are not computed at this stage, as they depend on the PLS
% results
%
% Inputs:
% - Ub_vect      : B x S x P matrix, S is #LCs, P is #bootstrap samples,
%                  bootstrapped behavior saliences for all LCs
% - Vb_vect      : M x S xP matrix, bootstrapped brain saliences for all LCs
% - boot_results : struct containing resluts from bootstrapping
%       - .Lxb,.Lyb,.LC_img_loadings_boot,.LC_behav_loadings_boot :
%           bootstrapping scores (see myPLS_get_PLSscores for details)
% - LC_behav_loadings_boot : B x S x P matrix, S is #LCs, 
%                            P is # bootstrap samples, bootstrapped behavior loadings
% - LC_img_loadings_boot   : M x S x P matrix, bootstrapped imaging loadings
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


%% Behavior and imaging saliences
% Computing mean, std and CIs for bootstrap saliences
boot_stats.Ub_mean=mean(Ub_vect,3); % mean
boot_stats.Ub_std=std(Ub_vect,[],3); % standard deviation
boot_stats.Ub_lB=prctile(Ub_vect,2.5,3); % lower bound of 95% CI
boot_stats.Ub_uB=prctile(Ub_vect,97.5,3); % lower bound of 95% CI

boot_stats.Vb_mean=mean(Vb_vect,3); % mean
boot_stats.Vb_std=std(Vb_vect,[],3); % standard deviation
boot_stats.Vb_lB=prctile(Vb_vect,2.5,3); % lower bound of 95% CI
boot_stats.Vb_uB=prctile(Vb_vect,97.5,3); % lower bound of 95% CI


%% Behavior and imaging loadings
% Computing mean, std and CIs for bootstrap saliences
boot_stats.LC_behav_loadings_mean=mean(boot_results.LC_behav_loadings_boot,3); % mean
boot_stats.LC_behav_loadings_std=std(boot_results.LC_behav_loadings_boot,[],3); % standard deviation
boot_stats.LC_behav_loadings_lB=prctile(boot_results.LC_behav_loadings_boot,2.5,3); % lower bound of 95% CI
boot_stats.LC_behav_loadings_uB=prctile(boot_results.LC_behav_loadings_boot,97.5,3); % lower bound of 95% CI

boot_stats.LC_img_loadings_mean=mean(boot_results.LC_img_loadings_boot,3); % mean
boot_stats.LC_img_loadings_std=std(boot_results.LC_img_loadings_boot,[],3); % standard deviation
boot_stats.LC_img_loadings_lB=prctile(boot_results.LC_img_loadings_boot,2.5,3); % lower bound of 95% CI
boot_stats.LC_img_loadings_uB=prctile(boot_results.LC_img_loadings_boot,97.5,3); % lower bound of 95% CI


% %%% Behavior loadings
% 
% % Compute mean, std deviation, z-scores, p-values of behavior loadings
% for iter_lc = 1:size(signif_LC,1)   
%     for iter_behav = 1:nBehav
% 
%         this_lc = signif_LC(iter_lc);
%         
%         % Std across samples
%         std_behav_boot(iter_behav,iter_lc) = std(LC_behav_loadings_boot(iter_behav,iter_lc,:));
%                
%         % Z-score (original loading / std loading across samples)
%         zscore_behav_boot(iter_behav,iter_lc) = LC_behav_loadings(iter_behav,this_lc) / std_behav_boot(iter_behav,iter_lc);
%         
%         % P-values for z-scores
%         if zscore_behav_boot(iter_behav,iter_lc) >= 0
%             pvals_behav_boot(iter_behav,iter_lc) = 1-cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
%         elseif zscore_behav_boot(iter_behav,iter_lc) < 0
%             pvals_behav_boot(iter_behav,iter_lc) = cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
%         end
%         
%     end    
% end
% 
% %%% Imaging loadings
% 
% % Compute mean, std deviation, z-scores, p-values of imaging loadings
% 
% for iter_lc = 1:size(signif_LC,1)
%     for iter_img = 1:nImaging
%         
%         % Std across samples
%         std_img_boot(iter_img,iter_lc) = std(LC_img_loadings_boot(iter_img,iter_lc,:));
%         
%         % Z-score (original loading / std loading across samples)
%         zscore_img_boot(iter_img,iter_lc) = LC_img_loadings(iter_img,iter_lc) / std_img_boot(iter_img,iter_lc);
%         
%         % P-values for z-scores
%         if zscore_img_boot(iter_img,iter_lc) >= 0
%             pvals_img_boot(iter_img,iter_lc) = 1-cdf('norm',zscore_img_boot(iter_img,iter_lc),0,1);
%         elseif zscore_img_boot(iter_img,iter_lc) < 0
%             pvals_img_boot(iter_img,iter_lc) = cdf('norm',zscore_img_boot(iter_img,iter_lc),0,1);
%         end
%     end
% end
