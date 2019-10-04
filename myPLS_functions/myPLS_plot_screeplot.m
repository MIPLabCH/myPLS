function myPLS_plot_screeplot(S,Sp_vect,pvals_LC,save_opts)

% This function plots singular values (observed & surrogates obtained 
% with permutation testing), as well as the covariance explained by each
% latent component (LC).
%
% Inputs:
% - S           : L x L matrix, L is #LCs, singular values (diagonal matrix)
% - Sp_vect     : L x #permutations matrix, permuted singular values         
% - pvals_LC    : L x 1 vector with p-values for each LC
% - save_opts   : Options for result saving and plotting
%       - .output_path : output directory where figures are saved
%       - .prefix      : prefix of all results files (optional)


nLC = size(diag(S),1);
file_name = fullfile(save_opts.output_path,[save_opts.prefix '_LC_pvals_explainedCovariance']);
    
% Get cumulative amount of covariance explained by each LC
cumul_explCov = cumsum(diag(S).^2 / sum(diag(S).^2)) * 100;

% Plot singular values (H1) & explained variance (H2)
figure;
[AX,H1,H2] = plotyy(1:nLC,diag(S),1:nLC,cumul_explCov,'plot','plot');

% Set options for plots
AX(2).YLim = [0 100];
set(H1,'Marker','o');
set(H1,'Color','k');
set(H2,'Color',[0 0.4470 0.7410]); 
set(AX,{'Ycolor'},{'k';[0 0.4470 0.7410]}) 
set(get(AX(1),'XLabel'),'String','Latent components','FontSize',14);
set(get(AX(1),'Ylabel'),'String','Singular values','FontSize',14) 
set(get(AX(2),'Ylabel'),'String','Explained covariance','FontSize',14) 
set(AX(1),'Ytick',0:1000:AX(1).YLim(2)); 
set(AX(2),'Ytick',0:20:100); 
set(gcf,'Color','w'); % Set background as white
hold on

% Add mean (SD) permuted S for each LC
errorbar(1:nLC, mean(Sp_vect,2), std(Sp_vect,[],2),'r-'); 
legend({'observed','surrogates'},'FontSize',12);

% Display p-values for each LC
for iLC = 1:nLC
    str{iLC} = sprintf('LC%d (p=%5.3f)',iLC,pvals_LC(iLC));
    fprintf('%s\n',str{iLC});
end

set(gca,'XTick',1:nLC,'XTickLabel',str,'Box','off','TickDir','out');

title('Explained covariance by each LC','FontSize',16);
    
% Save figure
saveas(gcf,[file_name '.jpg']);
