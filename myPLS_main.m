% This is the main script of to run the myPLS Toolbox
% 
% ------------------------------STEPS-------------------------------------
%
% The script includes: 
%   1. Call of a script with PLS inputs and their description
%   2. Call of the functions to run PLS and plot the results
%
%
% --------------------------REQUIREMENTS----------------------------------
%
%
%
%
% ------------------------- ON PLS METHOD---------------------------------
% For general descriptions of PLS for medical image analysis, we refer to:
%
% ﻿- Krishnan, A., Williams, L.J., McIntosh, A.R., Abdi, H., 2011. Partial 
%    Least Squares (PLS) methods for neuroimaging: A tutorial and review. 
%    Neuroimage 56, 455–475. doi:10.1016/j.neuroimage.2010.07.034
%
%  - McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of 
%    neuroimaging data: Applications and advances. Neuroimage 23, 250–263. 
%    doi:10.1016/j.neuroimage.2004.07.020
%
%
% ----------------------------EXAMPLES------------------------------------
% Example applications, which used this script:
%
%  - Behavior PLS with brain network-based measures:
%
%     * Zöller, D., Sandini, C., Karahanoğlu, F.I., Padula, M.C., Schaer, 
%       M., Eliez, S., Van De Ville, D., 2019. Large-scale brain network 
%       dynamics provide a measure of psychosis and anxiety in 22q11.2 
%       deletion syndrome. Biol. Psychiatry Cogn. Neurosci. Neuroimaging 
%       in press. doi:10.1016/j.bpsc.2019.04.004
%
%  - PLS for multivariate analysis of group differences and developmental effects:
%
%     * Zöller, D., Schaer, M., Scariati, E., Padula, M.C., Eliez, S., Van 
%       De Ville, D., 2017. Disentangling resting-state BOLD variability 
%       and PCC functional connectivity in 22q11.2 deletion syndrome. 
%       Neuroimage 149, 85–97. doi:10.1016/j.neuroimage.2017.01.064
%
%     * Zöller, D., Padula, M.C., Sandini, C., Schneider, M., Scariati, 
%       E., Van De Ville, D., Schaer, M., Eliez, S., 2018. Psychotic 
%       symptoms influence the development of anterior cingulate BOLD 
%       variability in 22q11.2 deletion syndrome. Schizophr. Res. 193, 
%       319–328. doi:10.1016/j.schres.2017.08.003
%
%
% -----------------------------CREDITS------------------------------------
% Code written by Prof. Dimitri Van De Ville, Daniela Zoeller and Valeria
% Kebets, with subfunctions borrowed from PLS toolbox by Rotman Baycrest
% (https://www.rotman-baycrest.on.ca/index.php?section=84)
% 
% These scripts and functions are based on MyPLS scripts previously 
% published at https://miplab.epfl.ch/index.php/software/PLS
%
% 
% Please cite the following papers when using this code:
% 
% Zöller, D., Schaer, M., Scariati, E., Padula, M.C., Eliez, S., Van De 
% Ville, D., 2017. Disentangling resting-state BOLD variability and PCC 
% functional connectivity in 22q11.2 deletion syndrome. Neuroimage 149, 
% 85–97. doi:10.1016/j.neuroimage.2017.01.064
% 
% McIntosh, A.R., Lobaugh, N.J., 2004. Partial least squares analysis of 
% neuroimaging data: Applications and advances. Neuroimage 23, 250–263. 
% doi:10.1016/j.neuroimage.2004.07.020

addpath('./myPLS_functions')
addpath('./RotmanBaycrest')
addpath('./misc')


%% defining all the inputs
myPLS_inputs

%% run PLS analysis, including permutation testing and bootstrapping
res = myPLS_analysis(input,pls_opts);

%% saving the results and plotting the outputs
% create output directory if necessary
if ~exist(save_opts.output_path);mkdir(save_opts.output_path);end





%% Contribution of original variables to LVs
% ??? What is this needed for?
% % Brain & behavior structure coefficients (Correlations imaging/behavior variables - brain/behavior scores)
% 
% clear myBrainStructCoeff myBehavStructCoeff
% 
% % Brain structure coefficients
% for iter_lv = 1:numSignifLVs
%     this_lv = mySignifLVs(iter_lv);
%     
%     for iter_img = 1:size(X,2)
%         clear tmpy tmpx r p
%         tmpx = X(:,iter_img);
%         tmpy = Lx(:,this_lv);
%         [r,p] = corrcoef(tmpx,tmpy.');
%         myBrainStructCoeff(iter_img,iter_lv) = r(1,2);
%     end
%     
% end
% 
% % Behavior structure coefficients
% for iter_lv = 1:numSignifLVs
%     this_lv = mySignifLVs(iter_lv);
% 
%     for iter_behav = 1:size(Y,2),
%         clear tmpy tmpx r p
%         tmpx = Y(:,iter_behav);
%         tmpy = Ly(:,this_lv);        
%         [r,p] = corrcoef(tmpx,tmpy.');
%         myBehavStructCoeff(iter_behav,iter_lv) = r(1,2);
%     end
% end


