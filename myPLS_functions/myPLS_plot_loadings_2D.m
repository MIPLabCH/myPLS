function myPLS_plot_loadings_2D(loadings,nRois,signif_LC,out_dir)
%
% This function plots 2-dimensional loadings (e.g. correlation matrix)
% of all significant latent components (LCs)
%
% Inputs:
% - loadings    : M x L matrix, M is #imaging values (when vectorized), L
%                 is #components
% - nRois       : number of regions of interest (ROIs) in correlation matrix
% - signif_LC   : significant latent components to plot
% - out_dir     : output directory where figures are saved
       
cd(out_dir);

for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    % Convert to symmetric correlation matrix
    CM = jVecToSymmetricMat(loadings(:,this_lc),nRois);
    
    figure; imagesc(CM);
    colormap('jet'); colorbar;
    xlabel('ROIs');
    ylabel('ROIs');
    set(gca,'TickLabelInterpreter','none','FontSize',6,'Box','off');
    set(gcf,'Color','w');
    title(['LC' num2str(this_lc) ' - FC loadings']);
    name_fig = ['LC' num2str(this_lc) '_FC_loadings.jpg'];
    saveas(gcf,name_fig);
    
end
