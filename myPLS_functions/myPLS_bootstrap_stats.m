function [std_behav_boot,zscore_behav_boot,pvals_behav_boot,std_img_boot,zscore_img_boot,pvals_img_boot] = ...
    myPLS_bootstrap_stats(LC_behav_loadings,LC_behav_loadings_boot,LC_img_loadings,LC_img_loadings_boot,signif_LC)
%
% This function computes the mean and std deviation across the imaging &
% behavior loadings obtained with bootstrap resampling (myPLS_bootstrap_loadings.m).
% Z-scores & p-values are then calculated by dividing the original loadings
% by their bootstrap-estimated standard deviation. 
%
% Inputs:
% - LC_behav_loadings      : B x L matrix, B is #behaviors, L is #latent
%                            components (LCs), behavior loadings
% - LC_behav_loadings_boot : B x S x P matrix, S is #significant LCs, 
%                            P is # bootstrap samples, bootstrapped behavior loadings 
% - LC_img_loadings        : M x L matrix, M is #imaging, imaging loadings 
% - LC_img_loadings_boot   : M x S x P matrix, bootstrapped imaging 
%                            loadings for significant LCs
% - signif_LC              : significant latent components to consider (e.g. [1,2])
%
% Outputs:
% std_behav_boot           : B x S matrix, standard deviation of behavior loadings
%                            across bootstrap samples for significant LCs
% zscore_behav_boot        : B x S matrix, zscore of behavior loadings
%                            across bootstrap samples for significant LCs
% pvals_behav_boot         : B x S matrix, p-values of behavior loadings
% std_img_boot             : M x S matrix, standard deviation of imaging
%                            loadings across bootstrap samples for significant LCs
% zscore_img_boot          : M x S matrix, zscore of imaging loadings
%                            across bootstrap samples for significant LCs
% pvals_img_boot           : 18 x 18 x S matrix, p-values of imaging loading 
%                            for significant LCs
%

nBehav = size(LC_behav_loadings,1);
nImaging = size(LC_img_loadings,1); 
nBootstraps = size(LC_behav_loadings_boot,3);

%%% Behavior loadings

% Compute mean, std deviation, z-scores, p-values of behavior loadings
for iter_lc = 1:size(signif_LC,1)   
    for iter_behav = 1:nBehav

        this_lc = signif_LC(iter_lc);
        
        % Std across samples
        std_behav_boot(iter_behav,iter_lc) = std(LC_behav_loadings_boot(iter_behav,iter_lc,:));
               
        % Z-score (original loading / std loading across samples)
        zscore_behav_boot(iter_behav,iter_lc) = LC_behav_loadings(iter_behav,this_lc) / std_behav_boot(iter_behav,iter_lc);
        
        % P-values for z-scores
        if zscore_behav_boot(iter_behav,iter_lc) >= 0
            pvals_behav_boot(iter_behav,iter_lc) = 1-cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        elseif zscore_behav_boot(iter_behav,iter_lc) < 0
            pvals_behav_boot(iter_behav,iter_lc) = cdf('norm',zscore_behav_boot(iter_behav,iter_lc),0,1);
        end
        
    end    
end

%%% Imaging loadings

% Compute mean, std deviation, z-scores, p-values of imaging loadings

for iter_lc = 1:size(signif_LC,1)
    for iter_img = 1:nImaging
        
        % Std across samples
        std_img_boot(iter_img,iter_lc) = std(LC_img_loadings_boot(iter_img,iter_lc,:));
        
        % Z-score (original loading / std loading across samples)
        zscore_img_boot(iter_img,iter_lc) = LC_img_loadings(iter_img,iter_lc) / std_img_boot(iter_img,iter_lc);
        
        % P-values for z-scores
        if zscore_img_boot(iter_img,iter_lc) >= 0
            pvals_img_boot(iter_img,iter_lc) = 1-cdf('norm',zscore_img_boot(iter_img,iter_lc),0,1);
        elseif zscore_img_boot(iter_img,iter_lc) < 0
            pvals_img_boot(iter_img,iter_lc) = cdf('norm',zscore_img_boot(iter_img,iter_lc),0,1);
        end
    end
end
