function myPLS_plot_results(res,save_opts)

% Function to create and save the output plots for PLS analysis
%
% INPUTS:
%   - res: Struct containing results from PLS - required fields for plots:
%     - .X0, .Y0          : unnormalized input matrices
%     - .X, .Y            : normalized input matrices
%     - .design_names     : names of design variables (have only changed in
%                           case of contrast PLS)
%     - .grouping,.group_names : copied from input for plotting
%     - .U                : B x L matrix, L is #latent components (LC), behavior saliences
%     - .V                : M x L matrix, imaging saliences
%     - .S                : L x L matrix, singular values (diagonal matrix)
%     - .explCovLC        : covariance explained by each LC
%     - .LC_pvals         : p-value for each LC (from permutation testing)
%     - .Lx               : brain scores
%     - .Ly               : behavior/design scores
%     - .LC_img_loadings  : corr(Lx,X)
%     - .LC_behav_loadings: corr(Ly,Y)
%     - .boot_results : struct with bootstrapping results
%           - .Ub_vect     : 3D matrix with bootstrapping samples of U
%           - .Vb_vect     : 3D matrix with bootstrapping samples of V
%           - .Lxb,.Lyb,.LC_img_loadings_boot,.LC_behav_loadings_boot :
%               3D matrices with bootstrapping PLS scores
%   - save_opts: Options for result saving and plotting
%       - .output_path   : path where to save the results
%       - .prefix        : prefix of all results files
%       - .img_type      : Specify how to plot the results
%              'volume' for voxel-based data in nifti Format - results 
%                       will be displayed as bootstrap ratios in a brain map
%              'corrMat' for ROI-to-ROI correlation matrix - results will 
%                       be displayed as bootstrap ratios in a correlation matrix
%              'barPlot' for any type of brain data in already vectorized 
%                       form - results will be displayed as barplots
%       - .mask_file     : gray matter mask, only required if imagingType='volume'
%       - .grouped_plots : binary variable indicating if groups should be 
%                          considered during plotting
%              0 = plotting ignoring grouping
%              1 = plotting cosidering grouping


% remove grouping information if it should be ignored for plotting
if ~save_opts.grouped_plots
    res.grouping=ones(size(res.grouping));
end

%% find significant latent variables
signif_LC=find(res.LC_pvals<save_opts.alpha);

%% plot the brain scores
disp('... Plotting brain scores ...')
myPLS_plot_subjScores(res.Lx,res.Ly,res.group_names,res.grouping,signif_LC);
print(gcf,fullfile(save_opts.output_path,[save_opts.prefix '_corrLxLy']),'-depsc2','-painters');
disp(' ');

%% plot imaging saliences and loadings
switch save_opts.img_type
    case 'volume'
        disp('... Saving bootstrap ratio maps ...')
        bootstrap_ratios=res.V./res.boot_results.Vb_std;
        myPLS_plot_loadings_3D(bootstrap_ratios,'BSR',signif_LC,...
            save_opts.BSR_thres(2),save_opts.BSR_thres(1),save_opts);
        disp(' ')
        
        disp('... Saving imaging loading maps ...')
        myPLS_plot_loadings_3D(res.LC_img_loadings,'img_loadings',signif_LC,...
            save_opts.load_thres(2),save_opts.load_thres(1),save_opts);
        disp(' ')
    case 'corrMat'
        disp('... Saving bootstrap ratio matrix ...')
        bootstrap_ratios=res.V./res.boot_results.Vb_std;
        myPLS_plot_loadings_2D(bootstrap_ratios,'BSR',signif_LC,...
            save_opts.BSR_thres(2),save_opts.BSR_thres(1),save_opts);
        disp(' ')
        
        disp('... Saving imaging loadings matrix ...')
        myPLS_plot_loadings_2D(res.LC_img_loadings,'img_loadings',signif_LC,...
            save_opts.load_thres(2),save_opts.load_thres(1),save_opts);
        disp(' ')
end

%% save results struct
disp('... Saving results ...')
save(fullfile(save_opts.output_path,[save_opts.prefix '_res']),'res','-v7.3');
disp(' ');
