function myPLS_plot_loadings_3D(loadings,mask_file,signif_LC,pos_threshold,neg_threshold,out_dir)
%
% This function displays 3-dimensional loadings (e.g. volume) of all 
% significant latent components (LCs). The loadings are constrained
% to a binary mask (e.g., DARTEL group template) and thresholded. 
% Writes loadings onto volume (.nii) for each significant LC.
% Requires SPM for reading/writing loadings onto volume, and Slover to
% display loadings.
%
% Inputs:
% - loadings       : M x L matrix, M is #imaging values (when vectorized),
%                     L is #components
% - mask_file      : binary mask filename (.nii or .img). Has to be in the
%                     same dimensions as loadings.
% - signif_LC      : significant latent components to plot
% - pos_threshold  : threshold for loadings' visualization (positive values)
% - neg_threshold  : threshold for loadings' visualization (negative values)
% - out_dir        : output directory where figures are saved

cd(out_dir);

% Load mask
mask = spm_read_vols(mask_file);
mask_idx = find(mask ~= 0);

% Make sure the mask is binary
a = unique(mask);
if size(a,1)>2,
    disp('Please input binary mask');
end
clear a
        
                
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    %%% 1. Write loadings onto volume
    
    % Constrain BSR to mask
    loadings_3D = zeros(size(mask));
    loadings_3D(mask_idx) = loadings(:,this_lc);
    Vi = mask_file;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = [out_dir '/LC' num2str(this_lc) '_img_loadings.nii'];
    spm_write_vol(Vi,loadings_3D);
    
    %%% 2. Display loadings
    
    % Load slover object
    load(fullfile(scripts_dir,'myobj_axial.mat')); % can be changed for 
    % sagittal/coronal view
    
    % Load loadings volume to get min/max values
    clear Ai A myMin myMax
    Ai = spm_vol([out_dir '/LC' num2str(this_lc) '_img_loadings.nii']);
    A = spm_read_vols(Ai);
    this_max = max(max(max(A)));
    this_min = min(min(min(A)));
    
    % Overlay loadings on template brain
    figure(11*iter_lc);
    set(gcf,'name',['LC' num2str(this_lc) ' - Imaging loadings']);
    set(gcf,'Position',[484 77 560 751]);
    myobj.img(1).vol = spm_vol(fullfile(scripts_dir,'ch2.nii')); % template
    myobj.img(2).vol = spm_vol([out_dir '/LC' num2str(this_lc) '_img_loadings.nii']); % loadings' positive values
    myobj.img(3).vol = spm_vol([out_dir '/LC' num2str(this_lc) '_img_loadings.nii']); % loadings' negative values
    myobj.img(2).range = [pos_threshold this_max];
    myobj.img(3).range = [neg_threshold this_min];
    myobj.figure = 11*iter_lc ;
    paint(myobj);
    
end