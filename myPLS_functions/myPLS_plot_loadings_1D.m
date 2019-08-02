function myPLS_plot_loadings_1D(loadings,std_loadings,var_type,var_names,signif_LC,nGroups_PLS,out_dir)
%
% This function plots one-dimensional loadings (e.g. behavior, imaging
% metrics) of all significant latent components (LCs) as barcharts.
%
% Inputs:
% - loadings       : B x L matrix, B is #variables, L is #components
% - std_loadings   : B x L matrix, standard deviation of loadings
% - var_type       : string, type of variable plotted
% - var_names      : string, names of varibles
% - signif_LC      : significant LCs to plot (e.g. [1,2])
% - nGroups_PLS    : number of groups used for PLS analysis, will determine
%                    if each group has its own set of saliences
% - out_dir        : output directory where figures are saved

nVars = size(loadings,1);

cd(out_dir);

% Plot loadings
for iter_lc = 1:length(signif_LC)
    this_lc = signif_LC(iter_lc);
    
    figure;
    %bar(loadings(:,this_lc));
    bar(reshape(loadings(:,this_lc),[nVars nGroups_PLS]));  
    errorbar((1:nVars)-.15+.30*(1-1), ...
    loadings(nVars*(1-1)+(1:CONST_NUM_BEHAV),this_lc),...
    std_loadings(nVars*(1-1)+(1:nVars),this_lc),'r.','MarkerSize',10);
    hold off
    xticks(1:nVars);
    xticklabels(var_names);
    set(gca,'TickLabelInterpreter','none','FontSize',6,'Box','off');
    set(gcf,'Color','w');
    xtickangle(45);
    xlabel([var_type ' variables']);
    ylabel('Correlation');
    title(['LC' num2str(this_lc) ' - ' var_type ' loadings']);
    name_fig = ['LC' num2str(this_lc) '_' var_type '_loadings.jpg'];
    saveas(gcf,name_fig);
end


