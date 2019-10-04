function myPLS_table_loadings(var_type,var_PLS_type,vars_PLS,vars_std,...
    vars_lB,vars_uB,save_opts,var_names,signif_LC)

% This function writes a .csv file containing a table reporting the 
% loadings & their standard error (standard deviation/confidence intervals)
%
% Inputs:
% - var_type       : string, type of variable plotted (e.g., 'Imaging','Design', 'Behavior')
% - var_PLS_type   : name of the PLS variable type for y-axis and for filenames
% - vars_PLS       : B x L matrix, B is #variables, L is #latent components
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

    switch save_opts.errorbar_mode
        case 'std'
            colNames = {'Loadings','SD'};
            A = [vars_PLS(:,iLC) vars_std(:,iLC)];

            % Uncomment to store table values as cell array to limit the number of decimals
            % + comment previous line
%             A = cell(1,1);
%             for iter_var = 1:numel(var_names)
%                 A{iter_var,1} = num2str(vars_PLS(iter_var,iLC),'%0.2f');
%                 A{iter_var,2} = num2str(vars_std(iter_var,iLC),'%0.2f');
%             end                        
            
        case 'CI'
            colNames = {'Loadings','Lower_CI','Upper_CI'};
            A = [vars_PLS(:,iLC) vars_lB(:,iLC) vars_uB(:,iLC)];
    end
        
    file_name = fullfile(save_opts.output_path,...
        [save_opts.prefix '_LC' num2str(this_lc) '_' var_type '_' var_PLS_type '.csv']);
    
    
    T = array2table(A,'VariableNames',colNames,'RowNames',rowNames);
    % Uncomment if table values were stored as cell array + comment previous line
    %T = cell2table(A,'VariableNames',colNames,'RowNames',rowNames);
    writetable(T,file_name);
end
