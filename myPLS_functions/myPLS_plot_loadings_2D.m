function myPLS_plot_loadings_2D(var,var_name,signif_LC,pos_threshold,neg_threshold,save_opts)

%
% This function plots 2-dimensional loadings (e.g. correlation matrix)
% of all significant latent components (LCs)
%
% Inputs:
% - var            : M x L matrix, M is #imaging values (when vectorized),
%                     L is #components
% - var_name       : name of the variable for filenames
% - signif_LC      : significant latent components to plot
% - pos_threshold  : threshold for loadings' visualization (positive values)
% - neg_threshold  : threshold for loadings' visualization (negative values)
% - save_opts      : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)

nImg=size(var,1);

% get number of regions (assumed that the vector contains vectorized
% elements of the upper triangular part of a symmetric matrix)
nRois=1;
while cumsum(1:nRois)<=nImg
    nRois=nRois+1;
end

for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    file_name = fullfile(save_opts.output_path,[save_opts.prefix '_LC' num2str(this_lc) '_' var_name]);
    
    % Convert to symmetric correlation matrix
    CM = jVecToSymmetricMat(var(:,this_lc),nRois,1);
    
    % thresholding
    CM(CM>0&CM<pos_threshold)=0;
    CM(CM<0&CM>neg_threshold)=0;
    
    % plotting
    figure; imagesc(CM);
    colormap('jet'); colorbar;
    xlabel('ROIs');
    ylabel('ROIs');
    set(gca,'TickLabelInterpreter','none','FontSize',6,'Box','off');
    set(gcf,'Color','w');
    title(['LC' num2str(this_lc) ' - ' var_name]);
    saveas(gcf,[file_name '.jpg']);
    
end
