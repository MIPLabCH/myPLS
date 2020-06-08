function myPLS_plot_results(res,save_opts)

% This function creates and saves the output plots of the PLS analysis
%
% Inputs:
% - res: Struct containing results from PLS - required fields for plots:
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
%           - .*_mean : mean of bootstrapping distributions
%           - .*_std : standard deviation of bootstrapping distributions
%           - .*_lB : lower bound of 95% confidence interval of bootstrapping distributions
%           - .*_uB : upper bound of 95% confidence interval of bootstrapping distributions
% - save_opts: Options for result saving and plotting
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
%       - .plot_boot_samples : binary variable indicating if bootstrap
%                          samples should be plotted in bar plots
%       - .errorbar_mode : 'std' = plotting standard deviations
%                          'CI' = plotting 95% confidence intervals
%       - .hl_stable	 : binary variable indicating if stable bootstrap
%                          scores should be highlighted 


% Remove grouping information if necessary for plotting
if ~save_opts.grouped_plots
    res.grouping = ones(size(res.grouping));
end

% Number and IDs of groups
groupIDs = unique(res.grouping);
nGroups = length(groupIDs);

% number of behavior scores
nBehav = size(res.Y,2);

%% Find significant latent components (LCs)

signif_LC = find(res.LC_pvals < save_opts.alpha);

%% Plot explained variance

myPLS_plot_screeplot(res.S,res.Sp_vect,res.LC_pvals,save_opts);

%% Plot null distribution of significant LCs

myPLS_plot_nulldistrib_Sp(res.S,res.Sp_vect,signif_LC,save_opts);

%% Plot PLS imaging & behavior/design scores

disp('... Plotting brain scores ...');

myPLS_plot_subjScores(res.Lx,res.Ly,res.group_names,res.grouping,signif_LC,save_opts);
disp(' ');

%% Plot PLS imaging saliences & loadings

switch save_opts.img_type
    case 'volume'
        disp('... Saving salience bootstrap ratio maps ...');
        bootstrap_ratios = res.V./res.boot_results.Vb_std;
        myPLS_plot_loadings_3D(bootstrap_ratios,'BSR_saliences',signif_LC,...
            save_opts.BSR_thres(2),save_opts.BSR_thres(1),save_opts);
        disp(' ')
        
        disp('... Saving loading bootstrap ratio maps ...');
        bootstrap_ratios = res.LC_img_loadings./res.boot_results.LC_img_loadings_std;
        myPLS_plot_loadings_3D(bootstrap_ratios,'BSR_loadings',signif_LC,...
            save_opts.BSR_thres(2),save_opts.BSR_thres(1),save_opts);
        disp(' ')
        
        disp('... Saving imaging loading maps ...');
        myPLS_plot_loadings_3D(res.LC_img_loadings,'img_loadings',signif_LC,...
            save_opts.load_thres(2),save_opts.load_thres(1),save_opts);
        disp(' ')
        
    case 'corrMat'
        disp('... Saving bootstrap ratio matrix ...');
        bootstrap_ratios = res.V./res.boot_results.Vb_std;
        myPLS_plot_loadings_2D(bootstrap_ratios,'BSR',signif_LC,...
            save_opts.BSR_thres(2),save_opts.BSR_thres(1),save_opts);
        disp(' ')
        
        disp('... Saving imaging loadings matrix ...');
        myPLS_plot_loadings_2D(res.LC_img_loadings,'img_loadings',signif_LC,...
            save_opts.load_thres(2),save_opts.load_thres(1),save_opts);
        disp(' ');
        
    case 'barPlot'
        disp('... Saving image salience bar plots ...');
        
        % Set up figure position
        if isfield(save_opts,'fig_pos_img') && ~isempty(save_opts.fig_pos_img)
            fig_pos = save_opts.fig_pos_img;
        else
            fig_pos = [440   378   560   420];
        end
        
        % Plot saliences
        % Set up vector for scatter plots
        if isfield(res.boot_results,'Vb_vect')
            vars_b_vect = res.boot_results.Vb_vect;
        else
            vars_b_vect = [];
        end
        % Call 1D bar plot function
        myPLS_plot_loadings_1D('Imaging','Saliences',res.V,vars_b_vect,...
            res.boot_results.Vb_mean,res.boot_results.Vb_std,...
            res.boot_results.Vb_lB,res.boot_results.Vb_uB,...
            res.img_names,signif_LC,1,fig_pos,save_opts);
        
        % Plot loadings
        % Set up vector for scatter plots
        if isfield(res.boot_results,'LC_img_loadings_boot')
            vars_b_vect = res.boot_results.LC_img_loadings_boot;
        else
            vars_b_vect = [];
        end
        % Call 1D bar plot function
        myPLS_plot_loadings_1D('Imaging','Loadings',res.LC_img_loadings,vars_b_vect,...
            res.boot_results.LC_img_loadings_mean,...
            res.boot_results.LC_img_loadings_std,...
            res.boot_results.LC_img_loadings_lB,...
            res.boot_results.LC_img_loadings_uB,...
            res.img_names,signif_LC,1,fig_pos,save_opts);
        disp(' ')
end


%% Plot behavior/design saliences & loadings (always as bar plots)

disp('... Saving image salience bar plots ...');
% Set up figure position
if isfield(save_opts,'fig_pos_behav') && ~isempty(save_opts.fig_pos_behav)
    fig_pos = save_opts.fig_pos_behav;
else
    fig_pos = [440   378   560   420];
end

% Set up variable type ('Behavior' or 'Design')
if ~isfield(save_opts,'behav_type') || isempty(save_opts.behav_type)
    var_type = 'Design'; % by default we call these variables 'Design'
else
    if strfind(save_opts.behav_type,'contrast')
        var_type = 'Design';
    else
        var_type='Behavior';
    end
end

% Set up variable names
var_names = repmat(res.design_names,nGroups,1);
if save_opts.grouped_plots == 1
    tmp = repmat(res.group_names,1,nBehav);
    for ii = 1:numel(var_names); var_names{ii} = [var_names{ii} ' (' tmp{ii} ')']; end;
    var_names = var_names(:);
end

% Plot saliences
% Set up vector for scatter plots
if isfield(res.boot_results,'Ub_vect')
    vars_b_vect = res.boot_results.Ub_vect;
else
    vars_b_vect = [];
end

% Call 1D bar plot function
myPLS_plot_loadings_1D(var_type,'Saliences',res.U,vars_b_vect,...
    res.boot_results.Ub_mean,res.boot_results.Ub_std,...
    res.boot_results.Ub_lB,res.boot_results.Ub_uB,...
    var_names,signif_LC,nGroups,fig_pos,save_opts);
              
% Plot loadings
% Set up vector for scatter plots
if isfield(res.boot_results,'LC_behav_loadings_boot')
    vars_b_vect = res.boot_results.LC_behav_loadings_boot;
else
    vars_b_vect = [];
end

% Call 1D bar plot function
myPLS_plot_loadings_1D(var_type,'Loadings',res.LC_behav_loadings,vars_b_vect,...
    res.boot_results.LC_behav_loadings_mean,...
    res.boot_results.LC_behav_loadings_std,...
    res.boot_results.LC_behav_loadings_lB,...
    res.boot_results.LC_behav_loadings_uB,...
    var_names,signif_LC,nGroups,fig_pos,save_opts);
disp(' ');
        
%% Save results struct

disp('... Saving results ...');
save(fullfile(save_opts.output_path,[save_opts.prefix '_res']),'res','-v7.3');
disp(' ');
