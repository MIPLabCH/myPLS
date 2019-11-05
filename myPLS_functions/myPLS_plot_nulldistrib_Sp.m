function myPLS_plot_nulldistrib_Sp(S,Sp_vect,signif_LC,save_opts)

% This function plots the null distribution of singular values obtained
% with permutation testing
%
% Inputs:
% - S          : L x L matrix, singular values (diagonal matrix)
% - Sp_vect    : matrix with permuted singular values
% - signif_LC      : significant latent components to plot
% - save_opts      : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)

S_vect=diag(S);

for iLC = 1:length(signif_LC)
    this_lc = signif_LC(iLC);
    
    this_S = S_vect(this_lc);
    this_Sp_vect = Sp_vect(this_lc,:);

    file_name = fullfile(save_opts.output_path,[save_opts.prefix '_LC' num2str(this_lc) '_nullDistribPermutedSingVals']);
    
    figure;
    histogram(this_Sp_vect,50,'FaceColor','b');
    title(['LC' num2str(this_lc) ' - Null distribution of singular values'],'FontSize',14);
    set(gca,'FontSize',12,'Box','off');
    set(gcf,'Color','w');
    a = std(this_Sp_vect);
    xlim([min(this_Sp_vect)-a this_S+a]);
    hold on
    line([this_S, this_S], ylim,'LineStyle',':','LineWidth',1.5,'Color', 'r');
    hold off
    
    saveas(gcf,[file_name '.jpg']);
 
end
