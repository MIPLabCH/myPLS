function [input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts)

% This function sets up the defaults and checks the validity of inputs for myPLS_analysis 

disp('... Initialization ...')

%% 1. Inputs

% Compatibility with X0 and Y0 inputs
if ~isfield(input,'brain_data') && isfield(input,'X0')
    input.brain_data = input.X0;
elseif ~isfield(input,'brain_data') && ~isfield(input,'X0')
    error('brain data input missing');
end

if ~isfield(input,'behav_data') && isfield(input,'Y0')
    input.behav_data = input.Y0;
elseif ~isfield(input,'behav_data') && ~isfield(input,'Y0')
    error('behavior data input missing');
end

% Check that dimensions of X & Y are correct
if(size(input.brain_data,1) ~= size(input.behav_data,1))
    error('Input arguments X and Y should have the same number of rows');
end


% Number and IDs of groups
groupIDs = unique(input.grouping);
nGroups = length(groupIDs);

% Number of behavior/design variables
nBehav = size(input.behav_data,2);

% Number of imaging variables
nImg = size(input.brain_data,2);

% Create defaults for group and behavior names, if not specified
if isfield(input,'group_names') && numel(input.group_names) < nGroups
    disp('!!! Fewer group names than groups - group names will be replaced by standard naming');
    input.group_names=[];
end
if ~isfield(input,'group_names') || isempty(input.group_names)
    input.group_names = cell(nGroups,1);
    for iG = 1:nGroups
        input.group_names{iG} = ['group ' num2str(groupIDs(iG))];
    end
end
if numel(input.group_names) > nGroups
    disp('!!! More group names than groups - only the first ones will be used:');
    input.group_names=input.group_names(1:nGroups);
    for iG=1:nGroups
        disp(['   group ' num2str(iG) ': "' input.group_names{iG} '"']);
    end
end

if isfield(input,'behav_names') && numel(input.behav_names) < nBehav
    disp('!!! Fewer behavior names than brahvior scores - behavior names will be replaced by standard naming');
    input.behav_names=[];
end
if ~isfield(input,'behav_names') || isempty(input.behav_names)
    input.behav_names = cell(nBehav,1);
    for iB = 1:nBehav
        input.behav_names{iB} = ['behavior ' num2str(iB)];
    end
end
if numel(input.behav_names) > nBehav
    disp('!!! More behavior names than behavior scores - only the first ones will be used:');
    input.behav_names=input.behav_names(1:nBehav);
    for iB=1:nBehav
        disp(['   behavior ' num2str(iB) ': "' input.behav_names{iB} '"']);
    end
end

if ~isfield(input,'behav_names') || isempty(input.behav_names)
    if isfield(save_opts,'img_type') && strcmp(save_opts.img_type,'barPlot')
        nImg = size(input.brain_data,2);
        for ii = 1:nImg
            input.img_names{ii,1} = ['img ' num2str(ii)]; 
        end
    end
end


%% 2. PLS options

% Default analysis type: behavior PLS
if ~isfield(pls_opts,'behav_type') || isempty(pls_opts.behav_type)
    pls_opts.behav_type = 'behavior';
elseif ~(strcmp(pls_opts.behav_type,'behavior') || ...
        strcmp(pls_opts.behav_type,'contrast') || ...
        strcmp(pls_opts.behav_type,'contrastBehav') || ...
        strcmp(pls_opts.behav_type,'contrastBehavInteract'))
    error('Invalid behav_type selected')
end

% Set up defaults for PLS grouping
if ~isfield(pls_opts,'grouped_PLS')
    disp('PLS grouping not specified, falling back to defaults:')
    if contains(pls_opts.behav_type,'contrast')
        pls_opts.grouped_PLS = 0;
        disp('   PLS includes contrast --> PLS without consideration of groups')
    else
        pls_opts.grouped_PLS = 1;
        disp('   PLS without contrast --> PLS under consideration of groups')
    end
end

% Check if normalization is compatible with contrast PLS
if contains(pls_opts.behav_type,'contrast')
    if pls_opts.normalization_behav == 2 || pls_opts.normalization_behav == 4 || ...
            pls_opts.normalization_img == 2 || pls_opts.normalization_img == 4
        error('Normalization within groups selected, but contrast in Y -> please change normalization to be done across all subjects!')
    end
    if pls_opts.grouped_PLS == 1
        error('Grouped PLS selected, but contrast in Y --> please change grouped_PLS variable to 0!')
    end
end

% Set up defaults for permutation/bootstrapping grouping
if ~isfield(pls_opts,'grouped_perm') || isempty(pls_opts.grouped_perm)
    disp('Permutations grouping not specified, falling back to defaults:')
    if contains(pls_opts.behav_type,'contrast')
        pls_opts.grouped_perm = 0;
        disp('   PLS includes contrast --> Permutations without consideration of groups')
    else
        pls_opts.grouped_perm = 1;
        disp('   PLS without contrast --> Permutations under consideration of groups')
    end
end

if ~isfield(pls_opts,'grouped_boot') || isempty(pls_opts.grouped_boot)
    disp('Permutations grouping not specified, falling back to defaults:')
    if contains(pls_opts.behav_type,'contrast')
        pls_opts.grouped_boot = 0;
        disp('   PLS includes contrast --> Bootstrapping without consideration of groups')
    else
        pls_opts.grouped_boot = 1;
        disp('   PLS without contrast --> Bootstrapping under consideration of groups')
    end
end

% Compatibility check for grouping in PLS/Permutations/Bootstrapping
if pls_opts.grouped_boot ~= pls_opts.grouped_perm
    disp('!!! Grouping option for permutations and bootstrapping not identical - make sure that this is what you would like to do!')
end
if pls_opts.grouped_PLS ~= pls_opts.grouped_perm
    disp('!!! Grouping option for PLS and permutations not identical - make sure that this is what you would like to do!')
end

% Set up defaults for Procrustes transform in bootstrapping
if ~isfield(pls_opts,'boot_procrustes_mod') || isempty(pls_opts.boot_procrustes_mod)
    disp('Procrustes transform mode for bootstrapping not specified, falling back to default:')
    pls_opts.boot_procrustes_mod = 1;
    disp('   Computation of Procrustes transform only using U (behavior saliences)')
end
if pls_opts.boot_procrustes_mod ~= 1 && pls_opts.boot_procrustes_mod ~= 2
    error('Invalid value in pls_opts.boot_procrustes_mod -> please set to either 1 or 2');
end

% Set up default for saving bootstrapping data  
if ~isfield(pls_opts,'save_boot_resampling') || isempty(pls_opts.save_boot_resampling)
    disp('No selection whether bootstrapping samples should be saved, falling back to default:')
    if nImg > 1000
        pls_opts.save_boot_resampling = 0;
        disp('   More than 1000 imaging dimensions -> no saving of bootstrap samples')
    else
        pls_opts.save_boot_resampling = 1;
    	disp('   Less than 1000 imaging dimensions -> bootstrap samples will be saved')
    end
end

if pls_opts.save_boot_resampling && nImg > 1000
    disp('!!! Bootstrap sample saving selected, but more than 1000 imaging dimensions - this can cause very large files!')
end

%% Saving & plotting options

% Set up default prefix (with some of the input parameters
if ~isfield(save_opts,'prefix') || isempty(save_opts.prefix)
    save_opts.prefix = sprintf('myPLS_TYPE%s_NORM%d%d',pls_opts.behav_type, ...
        pls_opts.normalization_img, pls_opts.normalization_behav);
end

% Check if an output directory has been specified
if ~isfield(save_opts,'output_path') || isempty(save_opts.output_path)
    error('Please specify an output directory to save results');
end

% Create output directory if necessary
if ~exist(save_opts.output_path,'dir')
    disp(['Creating output directory: ' save_opts.output_path])
    mkdir(save_opts.output_path);
end

% Set default alpha level
if ~isfield(save_opts,'alpha') || isempty(save_opts.alpha)
    save_opts.alpha = 0.05;
end

% Check imaging variable type
if ~isfield(save_opts,'img_type') || isempty(save_opts.img_type)
    error('Please specify the imaging type (volume/corrMat/barPlot)')
end

if strcmp(save_opts.img_type,'volume')
    if ~isfield(save_opts,'mask_file') || isempty(save_opts.mask_file) || ~exist(save_opts.mask_file,'file')
        error('Imaging type volume: specify a valid gray matter mask file!')
    end
    if ~isfield(save_opts,'struct_file') || isempty(save_opts.struct_file) || ~exist(save_opts.struct_file,'file')
        error('Imaging type volume: specify a valid structural file for image background!')
    end
    
    % check if mask dimensionality fits data
    % Load mask
    mask_hdr=spm_vol(save_opts.mask_file);
    mask = spm_read_vols(mask_hdr);
    mask_idx = logical(mask);
    nVoxMask=nnz(mask_idx);
    if nVoxMask ~= nImg
        error('!!! Invalid gray matter mask: different number of non-zero voxels in mask than in data! The mask should be identical to the mask used for masking brain data (X)')
    end
end

if strcmp(save_opts.img_type,'volume') || strcmp(save_opts.img_type,'corrMat')
    % Check for BSR threshold
    if ~isfield(save_opts,'BSR_thres') || isempty(save_opts.BSR_thres)
        save_opts.BSR_thres = [-3 3];
        disp('No BSR threshold specified, falling back to default threshold = +/-3')
    elseif length(save_opts.BSR_thres) == 1
        % Compatibility with single threshold input
        if save_opts.BSR_thres > 0
            save_opts.BSR_thres = [-save_opts.BSR_thres save_opts.BSR_thres];
        else 
            save_opts.BSR_thres = [save_opts.BSR_thres -save_opts.BSR_thres];
        end
    elseif length(save_opts.BSR_thres) > 2 || save_opts.BSR_thres(1) > 0 || save_opts.BSR_thres(2) < 0
        disp('!!! BSR threshold should have the following form: [negative_thres positive_thres], please check your inputs!');
    end
    
    % Check for loading threshold
    if ~isfield(save_opts,'load_thres') || isempty(save_opts.load_thres)
        save_opts.load_thres = [-0.1 0.1];
        disp('No loadings threshold specified, falling back to default threshold = +/-0.1')
    elseif length(save_opts.load_thres) == 1
        % Compatibility with single threshold input
        if save_opts.load_thres > 0
            save_opts.load_thres = [-save_opts.load_thres save_opts.load_thres];
        else 
            save_opts.load_thres = [save_opts.load_thres -save_opts.load_thres];
        end
    elseif length(save_opts.load_thres) > 2 || save_opts.load_thres(1) > 0 || save_opts.load_thres(2) < 0
        disp('!!! Loadings threshold should have the following form: [negative_thres positive_thres], please check your inputs!');
    end
end

% Check for compatibility with pls_opts
% Compatibility check for grouping in PLS and plotting
if pls_opts.grouped_PLS ~= save_opts.grouped_plots
    disp('!!! Grouping option for PLS and plotting not identical - make sure that this is what you would like to do!')
end


% Set boot samples plotting according to boot samples saving
if pls_opts.save_boot_resampling == 0
    if isfield(save_opts, 'plot_boot_samples') && save_opts.plot_boot_samples == 1
        disp('!!! Cannot plot boot samples if they are not saved, setting plot_boot_samples to 0!')
    end
    save_opts.plot_boot_samples = 0;
end

disp(' ')
