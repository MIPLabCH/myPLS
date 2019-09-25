function myPLS_table_loadings(var_type,var_type_name,vars_PLS,vars_std,...
    vars_lB,vars_uB,save_opts,var_names,signif_LC)

% This function writes a table reporting the loadings & their standard error
% (standard deviation/confidence intervals)
%
% Inputs:
% - var_type       : string, type of variable plotted (e.g., 'Imaging','Design', 'Behavior')
% - var_type_name  : name of the variable type for table filenames
% - vars_PLS       : B x L matrix, B is #variables, L is #components
% - vars_std       : B x L matrix, bootstrapping standard deviation
% - vars_lB        : B x L matrix, lower bound of 95% CIs
% - vars_uB        : B x L matrix, upper bound of 95% CIs
% - save_opts: Options for results saving and plotting
%       - .output_path   : path where to save the table
%       - .prefix        : prefix of all results files
%       - .errorbar_mode : 'std' = report standard deviations
%                          'CI'  = report 95% confidence intervals
% - var_names      : string, names of variables
% - signif_LC      : significant LCs to consider (e.g. [1,2])


rowNames = var_names;

for iLC = 1:length(signif_LC)
    this_lc = signif_LC(iLC);

    switch save_opts
        case 'std'
            colNames = {'Loadings','SD'};
            A = [vars_PLS(:,iLC) std_behav_boot(:,iLC)];
        case 'CI'
            colNames = {'Loadings','Lower CI','Upper CI'};
            A = [vars_PLS(:,iLC) vars_lB(:,iLC) vars_uB(:,iLC)];
    end
    
    file_name = fullfile(save_opts.output_path,...
        [save_opts.prefix '_LC' num2str(this_lc) '_' var_type '_' var_type_name '.xlsx']);
    
    T = array2table(A,'VariableNames',colNames,'RowNames',rowNames);
    writetable(T,file_name,'FileType','spreadsheet')
end
