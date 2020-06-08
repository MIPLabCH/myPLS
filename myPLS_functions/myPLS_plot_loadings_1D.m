function myPLS_plot_loadings_1D(var_type,var_PLS_type,vars_PLS,vars_b_vect,vars_mean,vars_std,...
    vars_lB,vars_uB,var_names,signif_LC,nGroups,fig_pos,save_opts)

% This function plots one-dimensional loadings (e.g. behavior, 1D imaging
% metrics) of all significant latent components (LCs) as barcharts.
%
% Inputs:
% - var_type       : string, type of variable plotted ('Imaging', 'Design' or 'Behavior')
% - var_PLS_type   : name of the PLS variable type for y-axis and for filenames
% - vars_PLS       : B x L matrix, B is #variables, L is #components
% - vars_b_vect    : B x L x P matrix, P is #bootstrap samples
% - vars_mean      : B x L matrix, bootstrapping mean
% - vars_std       : B x L matrix, bootstrapping standard deviation
% - vars_lB        : B x L matrix, lower bound of 95% CIs
% - vars_uB        : B x L matrix, upper bound of 95% CIs
% - var_names      : string, names of variables
% - signif_LC      : significant LCs to plot (e.g. [1,2])
% - nGroups        : number of groups used for PLS analysis, will determine
%                    if each group has its own set of saliences
% - fig_pos        : position of the figure to plot
% - save_opts: Options for results saving and plotting
%       - .output_path   : path where to save the results
%       - .prefix        : prefix of all results files
%       - .plot_boot_samples : binary variable indicating if bootstrap
%                          samples should be plotted in bar plots
%       - .errorbar_mode : 'std' = plotting standard deviations
%                          'CI' = plotting 95% confidence intervals
%       - .hl_stable	 : binary variable indicating if stable bootstrap
%                          scores should be highlighted
%       - .grouped_plots : binary variable indicating if groups should be 
%                          considered during plotting
%              0 = plotting ignoring grouping
%              1 = plotting cosidering grouping


nVars = size(vars_PLS,1);
nBootstraps = size(vars_b_vect,3);
nLCs = size(vars_PLS,2);

nBehav = nVars/nGroups;

% Set up plotting colors depending on var_type
if strcmp(var_type,'Imaging')
    plot_col = [0    0.5    0.5];
    if nGroups > 1
        error('Imaging data cannot be displayed for multiple groups');
    end
end

if strcmp(var_type,'Behavior') || strcmp(var_type,'Design')
    plot_col =  [0  0   0.5;
                 1  0.4 0.4;
                 0  0.5 0.5];
    plot_col = plot_col(1:nGroups,:);
    if nGroups > 3
        error('Multiple groups plotting only implemented up to 3 groups');
    end
end

% Re-order to group together the same behavior variables of different groups
iter = 1;
for iB = 1:nBehav
    re_order(iter:iter+nGroups-1) = iB:nBehav:nVars;
    iter = iter+nGroups;
end

% Plot loadings
for iLC = 1:length(signif_LC)
    this_lc = signif_LC(iLC);
    
    file_name = fullfile(save_opts.output_path,[save_opts.prefix '_LC' num2str(this_lc) '_' var_type '_' var_PLS_type]);
    
    figure('position',fig_pos);
    hold on;
    
    if save_opts.hl_stable
        % Check 'significance' and mark yellow the strong effeects
        max1 = max(abs(vars_lB(:,this_lc)));max2=max(abs(vars_uB(:,this_lc)));maxPatch=max(max1,max2)*1.1;
        
        for iC = 1:nVars % for every component
            if vars_mean(iC,this_lc) < 0
                if vars_uB(iC,this_lc) < 0 % significant
                    patch(iC+.5*[1,-1,-1,1],2*maxPatch*[-1,-1,1,1],[1 1 0.8],'edgecolor','none','LineStyle', 'none');
                end
            else
                if vars_lB(iC,this_lc) > 0 % significant
                    patch(iC+.5*[1,-1,-1,1],2*maxPatch*[-1,-1,1,1],[1 1 0.8],'edgecolor','none','LineStyle', 'none');
                end
            end
        end
        clearvars max*
        
        % draw 0 line on top (because the original x axis ended up below the patches)
        hl = refline(0,0);hl.Color='k';hl.LineWidth=1;
    end
    
    % Create barplot
    b = bar(vars_mean(re_order,this_lc),'FaceAlpha',.5); 
    b.FaceColor = 'flat';
    
    for iG = 1:nGroups
        b.CData(iG:nGroups:end,:) = repmat(plot_col(iG,:),size(b.CData,1)/nGroups,1);
    end
    
    % Create scatterplot, if asked for
    if save_opts.plot_boot_samples
        rand_x = [1:nVars]+0.08*randn(nBootstraps,1); % helper variable for scatterplots
        vars_tmp=squeeze(vars_b_vect(re_order,this_lc,:))'; % putting data in a vector
%         s=scatter(rand_x(:),vars_tmp(:),12,plot_col,'filled','MarkerFaceAlpha',.5);
        s=scatter(rand_x(:),vars_tmp(:),12,repmat(repelem(plot_col,nBootstraps,1),nBehav,1),...
                'filled','MarkerFaceAlpha',.5);
    end
    
    % Add errorbars
    switch save_opts.errorbar_mode
        case 'CI'
            h = ploterr([1:nVars],vars_mean(re_order,this_lc),[],...
                {vars_lB(re_order,this_lc),vars_uB(re_order,this_lc)},'k.', 'abshhxy', 0.2);
        case 'std'
            h = ploterr([1:nVars],vars_mean(re_order,this_lc),[],...
                vars_std(re_order,this_lc),'k.', 'abshhxy', 0.2);
    end
    set(h(1),'marker','none'); % remove marker
    set(h(2),'LineWidth',1.5);

    % Axes on top, set y axes limits, set x axes limits
    if exist('s','var')
        y_lim = [min(s.YData)-0.05 max(s.YData)+0.05];
    else
        y_dat = [h(1).YData h(2).YData];
        y_lim = [min(y_dat)-0.05 max(y_dat)+0.05];
    end
    xlims = get(gca,'xlim');
    set(gca,'Layer', 'Top','ylim',y_lim,'xlim',[xlims(1)+0.5 xlims(2)-0.5]); 
    
    % Labels & title
    xticks(1:nVars);
    xticklabels(var_names);
    set(gca,'TickLabelInterpreter','none','FontSize',6,'Box','off');
    set(gcf,'Color','w');
    xtickangle(45);
    xlabel([var_type ' variables']);
    ylabel(var_PLS_type);
    title(['LC' num2str(this_lc) ' - ' var_type ' ' var_PLS_type]);
    
    % Save figure
    saveas(gcf,[file_name '.jpg']);
    
    % Write loadings & standard error to table
    myPLS_table_loadings(var_type,var_PLS_type,vars_PLS,vars_std,...
    vars_lB,vars_uB,save_opts,var_names,signif_LC)    
        
end
