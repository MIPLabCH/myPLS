function [Y0,design_names,nDesignScores] = myPLS_getY(behavMode,...
    behaviorData,subj_grouping,groupnames,behav_names)

% number and IDs of groups
groupIDs=unique(subj_grouping);
nGroups=length(groupIDs);

% number of behavior scores
nBehav=size(behaviorData,2);

% number of subjects
nSubj = size(behaviorData,1); 


%% compute contrast, if necessary
% contrasts are implemented for up to 3 groups, for 3-group analysis see
% Zoller et al., Schizophrenia Research 2018
if contains(behavMode,'contrast')
    if nGroups==1
        error('There is only one group - contrasts are only implemented for 2 or 3 groups');
    elseif nGroups==2
        % grouping vectors
        group1ID=subj_grouping==groupIDs(1);nGroup1=nnz(group1ID);
        group2ID=subj_grouping==groupIDs(2);nGroup2=nnz(group2ID);

        % contrast
        c = zeros(nSubj,1);
        c(group1ID) = -1/nGroup1;
        c(group1ID) = 1/nGroup2;

        % othonormalize the contrast regressors
        c_orig=c;
        c = orth(c_orig);
        if sign(c(1))~=sign(c_orig(1)); c = -c; end % correct the sign if switched after norm 

        % specify the name of the contrast
        contrastName={['contrast (' groupnames{groupIDs(2)}...
            ' > ' groupnames{groupIDs(1)} ')']};
    elseif nGroups==3
        % grouping vectors
        group1ID=subj_grouping==groupIDs(1);nGroup1=nnz(group1ID);
        group2ID=subj_grouping==groupIDs(2);nGroup2=nnz(group2ID);
        group3ID=subj_grouping==groupIDs(2);nGroup3=nnz(group3ID);

        % contrast
        c = zeros(nSubj,2);
        c(group1ID,1) = -1/(nGroup1+nGroup2);
        c(group2ID,1) = -1/(nGroup1+nGroup2);
        c(group3ID,1) = 1/nGroup3;
        c(group1ID,2) = -1/nGroup1;
        c(group2ID,2) = 1/nGroup2;

        % othonormalize the contrast regressors
        c_orig=c;
        c(:,1) = orth(c_orig(:,1));
        c(:,2) = orth(c_orig(:,2));
        if sign(c(1,1))~=sign(c_orig(1,1)); c(:,1) = -c(:,1); end  % correct the sign if switched after norm 
        if sign(c(1,2))~=sign(c_orig(1,2)); c(:,2) = -c(:,2); end  % correct the sign if switched after norm 
        
        % specify the name of the contrast
        contrastName={['contrast (' groupnames{3} ...
            '>' groupnames{1}...
            '/' groupnames{2} ...
             ')'],...
            ['contrast (' groupnames{2}...
            '>' groupnames{1} ')']};
    else
        error('There are more than 3 groups - contrasts are only implemented for 2 or 3 groups');
    end
end


%% specify Y0 according to Mode


%% for behavior PLS, just keep the behavioral variables in Y0
switch behavMode
    case 'behavior'
        Y0=behaviorData;
        design_names=behav_names;
    case 'contrast'
        Y0=c;
        design_names=contrastName;
    case 'contrastBehav'
        % combine contrast and behavior
        Y0=[c,behaviorData];
        design_names = [contrastName, behav_names];
    case 'contrastBehavInteract'
        if nGroups==2
            % behavior and interaction
            Y0=[];
            for iBeh = 1:nBehav
                % behavior data
                Y0(:,2*iBeh-1)=behaviorData(iBeh);
                
                % normalization before interaction calculation
                Y0(:,2*iBeh)=zscore(Y0(:,2*iBeh-1)); 
                
                % contrast by behavior interaction
                Y0(:,2*iBeh)=Y0(:,2*iBeh).*c(:,1);
                
                % design names
                design_names{2*iBeh-1} = behav_names{iBeh};
                design_names{2*iBeh} = [behav_names{iBeh} '*' contrastName{1}];
            end
        elseif nGroups==3
            % behavior and interaction
            Y0=[];
            for iBeh = 1:length(behav_names)
                % behavior data
                Y0(:,3*iBeh-2)=behaviorData(iBeh);
                
                % normalization before interaction calculation
                Y0(:,3*iBeh-1)=-zscore(Y0(:,3*iBeh-2)); 
                Y0(:,3*iBeh)=zscore(Y0(:,3*iBeh-2));
                
                % contrast by behavior interaction
                Y0(:,3*iBeh-1)=Y0(:,3*iBeh-1).*c(:,1);
                Y0(:,3*iBeh)=Y0(:,3*iBeh).*c(:,2);
                
                % design names
                design_names{3*iBeh-2} = behav_names{iBeh};
                design_names{3*iBeh-1} = [behav_names{iBeh} '*' contrastName{1}];
                design_names{3*iBeh} = [behav_names{iBeh} '*' contrastName{2}];
            end
        end
        
        % combine contrast and behavior
        Y0=[c,Y0];
        design_names = [contrastName, design_names];
end

nDesignScores=size(Y0,2);


