function myPLS_plot_loadings_3D(var,var_name,signif_LC,pos_threshold,neg_threshold,save_opts)
%
% This function displays 3-dimensional loadings (e.g. volume) of all 
% significant latent components (LCs). The loadings are constrained
% to a binary mask (e.g., DARTEL group template) and thresholded. 
% Writes loadings onto volume (.nii) for each significant LC.
% Requires SPM for reading/writing loadings onto volume, and Slover to
% display loadings.
%
% Inputs:
% - var            : M x L matrix, M is #imaging values (when vectorized),
%                     L is #components
% - var_name       : name of the variable for filenames
% - signif_LC      : significant latent components to plot
% - pos_threshold  : threshold for loadings' visualization (positive values)
% - neg_threshold  : threshold for loadings' visualization (negative values)
% - save_opts      : Options for result saving and plotting
%       - .mask_file   : binary mask filename (.nii or .img). Has to be in
%                        the same dimensions as loadings.
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)
%       - .struct_file : filename of structural file for background volume
%                        to overlay the results on

% Load mask
mask_hdr=spm_vol(save_opts.mask_file);
mask = spm_read_vols(mask_hdr);
mask_idx = logical(mask);

% Make sure the mask is binary
a = unique(mask);
if size(a,1)>2,
    disp('Please input binary mask');
end
clear a
           
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    
    file_name = fullfile(save_opts.output_path,[save_opts.prefix '_LC' num2str(this_lc) '_' var_name]);
    
    %%% 1. Write loadings onto volume
    
    % Constrain BSR to mask
    var_3D = zeros(size(mask));
    var_3D(mask_idx) = var(:,this_lc);
    Vi = mask_hdr;
    Vi.dt = [spm_type('float32') 0];
    Vi.fname = [file_name '.nii'];
    spm_write_vol(Vi,var_3D);
    
    %%% 2. Display loadings
    
    % Load slover object
    switch save_opts.volume_orientation
        case 'axial'
            load('myobj_axial.mat');
            slices=-30:5:55;
        case 'sagittal'
            load('myobj_sagittal.mat');
            slices=-51:6:51;
        case 'coronal'
            load('myobj_coronal.mat');
            slices=-85:8:51;
    end
    
    % Load loadings volume to get min/max values
    clear Ai A myMin myMax
    Ai = spm_vol([file_name '.nii']);
    A = spm_read_vols(Ai);
    this_max = max(A(:));
    this_min = min(A(:));
    
    S = spm_read_vols(spm_vol(save_opts.struct_file));
    s_max=max(S(:));
    
    % Overlay loadings on template brain
    figure;
    set(gcf,'name',['LC' num2str(this_lc) ' - ' var_name]);
    set(gcf,'Position',[484 77 560 751]);
    myobj.img(1).vol = spm_vol(save_opts.struct_file); % template
    myobj.img(2).vol = spm_vol([file_name '.nii']); % loadings' positive values
    myobj.img(3).vol = spm_vol([file_name '.nii']); % loadings' negative values
    myobj.img(1).range=[0 0.7*s_max];
    myobj.img(2).range=[min(pos_threshold,this_max-0.1) this_max];
    myobj.img(3).range=[max(neg_threshold,this_min+0.1) this_min];
    myobj.figure = gcf;
    myobj.slices=slices;
    myobj.xslices=5;
    paint(myobj);
    print(gcf,file_name,'-depsc2','-painters');
    
end