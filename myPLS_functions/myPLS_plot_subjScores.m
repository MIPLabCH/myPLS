function corr_LxLy = myPLS_plot_subjScores(Lx,Ly,names_groups,grouping,signif_LC,save_opts)
%
% This function plots the imaging & behavior scores of all significant 
% latent components (LCs). Subjects' scores are shown according to their 
% diagnostic group. Correlations between imaging & behavior scores is displayed.
%
% Inputs:
% - Lx            : N x L matrix, N is #subjects, L#components, imaging scores
% - Ly            : N x L matrix, behavioral  scores
% - names_groups  : string, names of diagnostic groups
% - grouping      : N x 1 vector, subject (diagnostic) grouping 
%                   e.g. [1,1,2] = subjects 1&2 belong to group 1,
%                   subject 3 belongs to group 2.
% - signif_LC     : significant latent components to plot, e.g. [1,2]
%
% Outputs:
% - corr_LxLy     : correlation between imaging and behavior scores

% number and IDs of groups
groupIDs=unique(grouping);
nGroups=length(groupIDs);

colors = {'b','r','c','g','m','y','w','k'}; % Matlab colors
plot_colors = colors(1:nGroups); % select as many colors as groups

disp('Correlations between imaging and behavioral scores');
for iter_lc = 1:size(signif_LC,1)
    this_lc = signif_LC(iter_lc);
    
    figure('position',[440   541   327   257]);
    for iG = 1:nGroups
        plot(Lx(grouping==groupIDs(iG),this_lc),...
            Ly(grouping==groupIDs(iG),this_lc),[plot_colors{iG} '.'],'MarkerSize',20);
        hold on
    end
    hold off
    title(['LC' num2str(this_lc) ' - Correlations between imaging and behavioral scores']);
    
    
    legend(names_groups,'Location','southeast');
    xlabel('Imaging scores');
    ylabel('Behavior/Design scores');
    
    set(gcf,'Color','w');
    set(gca,'Box','off');
    
    corr_LxLy(iter_lc) = corr(Lx(:,this_lc),Ly(:,this_lc));
    disp(['LC' num2str(this_lc) ': r = ' num2str(corr_LxLy(iter_lc),'%0.2f')]);
    
    print(gcf,fullfile(save_opts.output_path,[save_opts.prefix '_LC' num2str(iter_lc) '_corrLxLy']),'-depsc2','-painters');
end