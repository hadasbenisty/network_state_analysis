function [outputMat_grabs_regular, outputMat_grabs_ridge] = ...
    returnCVRSquared(animal,k,dataFolderPath,imagingdatapath,spike2path,airpuff_trial,howManyFacePCs)
[parcels_names, parcels_region_labels, finalindex, region_lut] = get_allen_meta_parcels; %#ok<ASGLU>

%dbstop if warning
tic
%**************************************************************************
% airpuff_trial input is 1 if it is an airpuff trial, 0 if not
% howManyFacePCs: number of face PCs to include (any number 1-100)
%
% extract:
% 
% pupil trace
% face PC1
% wheel trace
% EEG (30-100 Hz)
%
% regress parcel-wise using regular regression: regress(y,X)
% compute cross-validated R squared using k=10 folds
%
% regress parcel-wise using ridge regression: ridgeMML(y,R_mat,'true')
% compute cross-validated R squared using k=10 folds
%
% outputs: matrix of parcel-wise Rsquared values for grabs and rcamp
%
% row 1 = all behavioral vars
% row 2 = pupil
% row 3 = face PC 1
% row 4 = wheel
% row 5 = EEG bandpower 30-100 Hz
%**************************************************************************

% *************************************************************************
% load in data
% *************************************************************************

%cd(dataFolderPath); % change directory
% detrending params
params.deterend.filtLen = 150;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';

    

if isfile(fullfile(imagingdatapath,'Ca_traces_spt_patch14_Allen_dfff.mat'))&&isfile(fullfile(dataFolderPath,'smrx_signals_v3.mat'))
    
dFoF_parcells = load(fullfile(imagingdatapath,'Ca_traces_spt_patch14_Allen_dfff.mat')); %<- parcels data
smrx_sigs = load(fullfile(dataFolderPath,'smrx_signals_v3.mat'));
% get those parcels from one hemi
blue_parcels=dFoF_parcells.parcels_time_trace;
%add dfof into the same folders later 
%smrx_sigs.timestamps.timaging=smrx_sigs.timestamps.timaging(params.deterend.filtLen/2:end);

%cut off time at the end if not aligned
if size(blue_parcels,2)==length(smrx_sigs.timestamps.timaging)
elseif size(blue_parcels,2)<length(smrx_sigs.timestamps.timaging)
    smrx_sigs.timestamps.timaging=smrx_sigs.timestamps.timaging(1:size(blue_parcels,2));
elseif size(blue_parcels,2)>length(smrx_sigs.timestamps.timaging)
    blue_parcels=blue_parcels(:,1:length(smrx_sigs.timestamps.timaging));
end

%equal length sessions input jan 26 2021
if size(blue_parcels,2)>1800&&length(smrx_sigs.timestamps.timaging)>1800
    blue_parcels=blue_parcels(:,1:1800);
    smrx_sigs.timestamps.timaging=smrx_sigs.timestamps.timaging(1:1800);
else
end
proc_filelist = (dir(fullfile(spike2path, '*_proc.mat')));
ProcFileName=fullfile(spike2path,proc_filelist.name);
proc_output = load(fullfile(ProcFileName));

% 1) pupil*****************************************************************

pupilfile = dir(fullfile(dataFolderPath, 'pupil_clean.mat'));
if  isempty(pupilfile)
    disp('no pupil file');
    return;
end

%% for each mouse, load spont and airpuff folders and perform correlations

if ~isfield(smrx_sigs.timing, 'airpuffstart')
    if isfield(smrx_sigs.timing, 'stimstart') % if this is a vis session, save vis as airpuff so we'll know to avoid these segmenents
        smrx_sigs.timing.airpuffstart = smrx_sigs.timing.stimstart;
        smrx_sigs.timing.airpuffend = smrx_sigs.timing.stimend;
    else
        smrx_sigs.timing.airpuffstart = [];
        smrx_sigs.timing.airpuffend = [];
    end
end
dat=load(fullfile(pupilfile.folder, pupilfile.name));
if ~isfield(dat, 'areaii')
    warning('No Pupil on face map');
    return;
end
pupil_Norm = dat.areaii;
if size(pupil_Norm,1)==length(smrx_sigs.timing.pupilcamstart(1:end))
elseif size(pupil_Norm,1)<length(smrx_sigs.timing.pupilcamstart(1:end))
    % if pupil video is smaller than cam ticks then subset cam ticks
    smrx_sigs.timing.pupilcamstart=smrx_sigs.timing.pupilcamstart(1:size(pupil_Norm,1));
    smrx_sigs.timing.pupilcamend=smrx_sigs.timing.pupilcamend(1:size(pupil_Norm,1));
elseif size(pupil_Norm,1)>length(smrx_sigs.timing.pupilcamstart(1:end))
    %if cam ticks are more subset pupil and facemap video
    pupil_Norm=pupil_Norm(1:length(smrx_sigs.timing.pupilcamstart),:);
    proc_output.proc.motSVD{1,1}=proc_output.proc.motSVD{1,1}(1:length(smrx_sigs.timing.pupilcamstart),:);
end
pupil_time=smrx_sigs.timing.pupilcamstart;
%remove NaNs

%find the smallest of facemap, pupil, and imaging
if pupil_time(1) > smrx_sigs.timestamps.timaging(1)
   % if imaging starts earlier than pupil
   toremove_img=find(smrx_sigs.timestamps.timaging<pupil_time(1));
   smrx_sigs.timestamps.timaging(toremove_img)=[];
   blue_parcels(:,toremove_img)=[];
elseif pupil_time(1) < smrx_sigs.timestamps.timaging(1)
    %interp will take care of this case
else   
end

pupil_interp = interp1(pupil_time,pupil_Norm,smrx_sigs.timestamps.timaging);
pupil_interp(isnan(pupil_interp))=0;
pupil_sig=zscore(pupil_interp);
% what lohani and moberly did
% pupil=proc_output.proc.pupil.area; 
% pupil_time = smrx_sigs.timing.pupilcamstart(1:length(pupil));
% % interpolate and z score
% pupil_interp = interp1(pupil_time,pupil,smrx_sigs.timestamps.timaging);
% pupil_sig = zscore(pupil_interp);

% 2) face PC1**************************************************************
wholeFaceSVD_time = pupil_time;

wholeFaceSVD_interp = zeros(length(smrx_sigs.timestamps.timaging), 100);
    for i = 1:100
        tmp = proc_output.proc.motSVD{1,1}(:,i);
        
        tmp_interp = interp1(wholeFaceSVD_time,tmp,smrx_sigs.timestamps.timaging);
        tmp_interp(isnan(tmp_interp))=0;
        wholeFaceSVD_interp(:,i) = zscore(tmp_interp);
    end
 face_PC1 = wholeFaceSVD_interp(:,1:howManyFacePCs);
    
%3 wheel speed cont********************************************************
wheel_speed = smrx_sigs.channels_data.wheelspeed;
wheel_time = (1:length(wheel_speed))/5000;
%interpolate
wheel_speed_interp = interp1(wheel_time,wheel_speed,smrx_sigs.timestamps.timaging);  
wheel = zscore(wheel_speed_interp);

%**************************************************************************
% remove normally nan parcels
blue_parcels([21:26 53:56],:) = []; % remove nan values
blue_parcels_z = zeros(size(blue_parcels,1), size(blue_parcels,2));

% z score
    for i = 1:size(blue_parcels,1)
        tmp = blue_parcels(i,:)-nanmean(blue_parcels(i,:))./nanstd(blue_parcels(i,:));
        blue_parcels_z(i,:) = tmp;
    end  

%**************************************************************************
% if an airpuff trial-> remove data around airpuffs
%**************************************************************************
if airpuff_trial
    % TIMESTAMPS OF AIRPUFFS (seconds) smrx_sigs.timing.allairpuffstart
    % remove 10 second before and 60 seconds after airpuff (can change)
    time_to_remove_before = 1*10; %s * 10 Hz
    time_to_remove_after = 6*10; %s * 10 Hz
    %lav changed to 1 and 6.
    
    %1 find indicies
    %2 make binary matrix indicating what to delete
    airpuff_delete = zeros(1,length(smrx_sigs.timestamps.timaging));

    for i = 1:length(smrx_sigs.timing.allairpuffstart)
    
        [ ~, ix ] = min( abs( smrx_sigs.timestamps.timaging-smrx_sigs.timing.allairpuffstart(i) ) );
        if (ix-time_to_remove_before)>=1&&(ix+time_to_remove_after)<=length(smrx_sigs.timestamps.timaging)
        airpuff_delete(ix-time_to_remove_before:ix+time_to_remove_after) = 1;
        else
        end
    end

    % if airpuff_delete vector is too long, cut off last entries:
    if length(airpuff_delete) > length(smrx_sigs.timestamps.timaging)
        airpuff_delete = airpuff_delete(1:length(smrx_sigs.timestamps.timaging));
    end

    %3 remove these inds from each behavioral variable:
    pupil_sig(airpuff_delete == 1) = [];

    face_PC1(airpuff_delete == 1,:) = [];

    wheel(airpuff_delete == 1) = [];

    blue_parcels_z(:,airpuff_delete == 1) = [];
    
    blue_parcels(:,airpuff_delete == 1) = [];
      
end    
    
%remove NaNs from detrending filter and realign
% blue_parcels = blue_parcels(:,params.deterend.filtLen/2:end); 
% blue_parcels_z = blue_parcels_z(:,params.deterend.filtLen/2:end); 
% wheel=wheel(params.deterend.filtLen/2:end);
% face_PC1=face_PC1(params.deterend.filtLen/2:end,:);
% pupil_sig=pupil_sig(params.deterend.filtLen/2:end);

%find periods when all parcels have NaNs, which is due to the detrending
%filter or sometimes 1-3 frames after. controls for removal of timestamps
%earlier possibly  due to alignment with pupil (thus the detrending filter
%may no longer apply)

toremove=all(~isnan(blue_parcels))==0; % find indices where there are missing NaNs for all parcels
%the indices are always the same as find(all(isnan(blue_parcels))==1) since
%nans are usually missing from all parcels at once

pupil_sig(toremove== 1) = [];

face_PC1(toremove== 1,:) = [];
wheel(toremove== 1)= [];

blue_parcels_z(:,toremove == 1) = [];

blue_parcels(:,toremove == 1) = [];



%**************************************************************************
% regression using y as parcel data and X as behavioral variables:
% 10 fold cross validation using cvpartition function
%**************************************************************************

%1) using all behavioral variables GRABS:

% regular regression:

rsquared_grabs_allVars = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
for parcel_num = 1:size(blue_parcels,1) % for each parcel
    y = blue_parcels_z(parcel_num,:)';
    X = [pupil_sig face_PC1 wheel];
    K = 10; % number of folds
    cv = cvpartition(size(y,1),'kfold',K);
    cv_rsquared_values = zeros(K,size(y,2));
    for i = 1:K
        % training and testing indices for this fold
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        % return beta dims on training set:
        
        betas = regress(y(trainIdx,:), X(trainIdx,:)); %get beta weights.
        % predict using testing data
        Vm = (X(testIdx,:) * betas); % reconstruct data
        all_rsquareds = zeros(1,size(y,2));
        for j = 1:size(y,2)
            predicted = Vm(:,j);
            actual = y(testIdx,:);
            actual = actual(:,j);
            % caluculate R squared
            SSR = sum((predicted-nanmean(actual)).^2);
            SSTotal = sum((actual-nanmean(actual)).^2);
            all_rsquareds(j) = SSR/SSTotal;
        end
        cv_rsquared_values(i,:) = all_rsquareds;
    end
    cvrsquared = mean(cv_rsquared_values);
    rsquared_grabs_allVars(parcel_num) = cvrsquared;
end

% Ridge regression:
% rsquared_grabs_allVars_ridge = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
% for parcel_num = 1:size(blue_parcels,1) % for each parcel
%     y = blue_parcels_z(parcel_num,:)';
%     X = [pupil_sig face_PC1 wheel];
%     K = 10; % number of folds
%     cv = cvpartition(size(y,1),'kfold',K);
%     cv_rsquared_values = zeros(K,size(y,2));
%     for i = 1:K
%         % training and testing indices for this fold
%         trainIdx = cv.training(i);
%         testIdx = cv.test(i);
%         % return beta dims on training set:
%         [~, betas] = ridgeMML(y(trainIdx,:), X(trainIdx,:), true); %get ridge penalties and beta weights.
%         % predict using testing data
%         Vm = (X(testIdx,:) * betas); % reconstruct data
%         all_rsquareds = zeros(1,size(y,2));
%         for j = 1:size(y,2)
%             predicted = Vm(:,j);
%             actual = y(testIdx,:);
%             actual = actual(:,j);
%             % caluculate R squared
%             SSR = sum((predicted-mean(actual)).^2);
%             SSTotal = sum((actual-mean(actual)).^2);
%             all_rsquareds(j) = SSR/SSTotal;
%         end
%         cv_rsquared_values(i,:) = all_rsquareds;
%     end
%     cvrsquared = mean(cv_rsquared_values);
%     rsquared_grabs_allVars_ridge(parcel_num) = cvrsquared;
% end

%**************************************************************************
%3) using pupil only~ GRABS:

%regular regression:
rsquared_grabs_pupil = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
for parcel_num = 1:size(blue_parcels,1) % for each parcel
    y = blue_parcels_z(parcel_num,:)';
    X = pupil_sig;
    K = 10; % number of folds
    cv = cvpartition(size(y,1),'kfold',K);
    cv_rsquared_values = zeros(K,size(y,2));
    for i = 1:K
        % training and testing indices for this fold
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        % return beta dims on training set:
        betas = regress(y(trainIdx,:), X(trainIdx,:)); %get beta weights.
        % predict using testing data
        Vm = (X(testIdx,:) * betas); % reconstruct data
        all_rsquareds = zeros(1,size(y,2));
        for j = 1:size(y,2)
            predicted = Vm(:,j);
            actual = y(testIdx,:);
            actual = actual(:,j);
            % caluculate R squared
            SSR = sum((predicted-mean(actual)).^2);
            SSTotal = sum((actual-mean(actual)).^2);
            all_rsquareds(j) = SSR/SSTotal;
        end
        cv_rsquared_values(i,:) = all_rsquareds;
    end
    cvrsquared = mean(cv_rsquared_values);
    rsquared_grabs_pupil(parcel_num) = cvrsquared;
end

% Ridge regression:
% rsquared_grabs_pupil_ridge = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
% for parcel_num = 1:size(blue_parcels,1) % for each parcel
%     y = blue_parcels_z(parcel_num,:)';
%     X = pupil_sig;
%     K = 10; % number of folds
%     cv = cvpartition(size(y,1),'kfold',K);
%     cv_rsquared_values = zeros(K,size(y,2));
%     for i = 1:K
%         % training and testing indices for this fold
%         trainIdx = cv.training(i);
%         testIdx = cv.test(i);
%         % return beta dims on training set:
%         [~, betas] = ridgeMML(y(trainIdx,:), X(trainIdx,:), true); %get ridge penalties and beta weights.
%         % predict using testing data
%         Vm = (X(testIdx,:) * betas); % reconstruct data
%         all_rsquareds = zeros(1,size(y,2));
%         for j = 1:size(y,2)
%             predicted = Vm(:,j);
%             actual = y(testIdx,:);
%             actual = actual(:,j);
%             % caluculate R squared
%             SSR = sum((predicted-mean(actual)).^2);
%             SSTotal = sum((actual-mean(actual)).^2);
%             all_rsquareds(j) = SSR/SSTotal;
%         end
%         cv_rsquared_values(i,:) = all_rsquareds;
%     end
%     cvrsquared = mean(cv_rsquared_values);
%     rsquared_grabs_pupil_ridge(parcel_num) = cvrsquared;
% end

%**************************************************************************

%**************************************************************************
%5) using Face only~ GRABS:
rsquared_grabs_face = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
for parcel_num = 1:size(blue_parcels,1) % for each parcel
    y = blue_parcels_z(parcel_num,:)';
    X = face_PC1;
    K = 10; % number of folds
    cv = cvpartition(size(y,1),'kfold',K);
    cv_rsquared_values = zeros(K,size(y,2));
    for i = 1:K
        % training and testing indices for this fold
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        % return beta dims on training set:
        betas = regress(y(trainIdx,:), X(trainIdx,:)); %get beta weights.
        % predict using testing data
        Vm = (X(testIdx,:) * betas); % reconstruct data
        all_rsquareds = zeros(1,size(y,2));
        for j = 1:size(y,2)
            predicted = Vm(:,j);
            actual = y(testIdx,:);
            actual = actual(:,j);
            % caluculate R squared
            SSR = sum((predicted-mean(actual)).^2);
            SSTotal = sum((actual-mean(actual)).^2);
            all_rsquareds(j) = SSR/SSTotal;
        end
        cv_rsquared_values(i,:) = all_rsquareds;
    end
    cvrsquared = mean(cv_rsquared_values);
    rsquared_grabs_face(parcel_num) = cvrsquared;
end

%Ridge regression:
% rsquared_grabs_face_ridge = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
% for parcel_num = 1:size(blue_parcels,1) % for each parcel
%     y = blue_parcels_z(parcel_num,:)';
%     X = face_PC1;
%     K = 10; % number of folds
%     cv = cvpartition(size(y,1),'kfold',K);
%     cv_rsquared_values = zeros(K,size(y,2));
%     for i = 1:K
%         %training and testing indices for this fold
%         trainIdx = cv.training(i);
%         testIdx = cv.test(i);
%         %return beta dims on training set:
%         [~, betas] = ridgeMML(y(trainIdx,:), X(trainIdx,:), true); %get ridge penalties and beta weights.
%         %predict using testing data
%         Vm = (X(testIdx,:) * betas); % reconstruct data
%         all_rsquareds = zeros(1,size(y,2));
%         for j = 1:size(y,2)
%             predicted = Vm(:,j);
%             actual = y(testIdx,:);
%             actual = actual(:,j);
%             %caluculate R squared
%             SSR = sum((predicted-mean(actual)).^2);
%             SSTotal = sum((actual-mean(actual)).^2);
%             all_rsquareds(j) = SSR/SSTotal;
%         end
%         cv_rsquared_values(i,:) = all_rsquareds;
%     end
%     cvrsquared = mean(cv_rsquared_values);
%     rsquared_grabs_face_ridge(parcel_num) = cvrsquared;
% end
% %}
%**************************************************************************
%**************************************************************************
%7) using wheel only~ GRABS:
rsquared_grabs_wheel = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
for parcel_num = 1:size(blue_parcels,1) % for each parcel
    y = blue_parcels_z(parcel_num,:)';
    X = wheel;
    K = 10; % number of folds
    cv = cvpartition(size(y,1),'kfold',K);
    cv_rsquared_values = zeros(K,size(y,2));
    for i = 1:K
        % training and testing indices for this fold
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
        % return beta dims on training set:
        betas = regress(y(trainIdx,:), X(trainIdx,:)); %get beta weights.
        % predict using testing data
        Vm = (X(testIdx,:) * betas); % reconstruct data
        all_rsquareds = zeros(1,size(y,2));
        for j = 1:size(y,2)
            predicted = Vm(:,j);
            actual = y(testIdx,:);
            actual = actual(:,j);
            % caluculate R squared
            SSR = sum((predicted-mean(actual)).^2);
            SSTotal = sum((actual-mean(actual)).^2);
            all_rsquareds(j) = SSR/SSTotal;
        end
        cv_rsquared_values(i,:) = all_rsquareds;
    end
    cvrsquared = mean(cv_rsquared_values);
    rsquared_grabs_wheel(parcel_num) = cvrsquared;
end

% % Ridge regression:
% rsquared_grabs_wheel_ridge = zeros(1,size(blue_parcels,1)); % initialize variable to hold parcel rsqrd vals
% for parcel_num = 1:size(blue_parcels,1) % for each parcel
%     y = blue_parcels_z(parcel_num,:)';
%     X = wheel;
%     K = 10; % number of folds
%     cv = cvpartition(size(y,1),'kfold',K);
%     cv_rsquared_values = zeros(K,size(y,2));
%     for i = 1:K
%         % training and testing indices for this fold
%         trainIdx = cv.training(i);
%         testIdx = cv.test(i);
%         % return beta dims on training set:
%         [~, betas] = ridgeMML(y(trainIdx,:), X(trainIdx,:), true); %get ridge penalties and beta weights.
%         % predict using testing data
%         Vm = (X(testIdx,:) * betas); % reconstruct data
%         all_rsquareds = zeros(1,size(y,2));
%         for j = 1:size(y,2)
%             predicted = Vm(:,j);
%             actual = y(testIdx,:);
%             actual = actual(:,j);
%             % caluculate R squared
%             SSR = sum((predicted-mean(actual)).^2);
%             SSTotal = sum((actual-mean(actual)).^2);
%             all_rsquareds(j) = SSR/SSTotal;
%         end
%         cv_rsquared_values(i,:) = all_rsquareds;
%     end
%     cvrsquared = mean(cv_rsquared_values);
%     rsquared_grabs_wheel_ridge(parcel_num) = cvrsquared;
% end

%**************************************************************************

%**************************************************************************
% combine into one regular regression output matrix
outputMat_grabs_regular = vertcat(rsquared_grabs_allVars,rsquared_grabs_pupil,rsquared_grabs_face,rsquared_grabs_wheel);    

% combine into one ridge regression output matrix 
%outputMat_grabs_ridge = vertcat(rsquared_grabs_allVars_ridge,rsquared_grabs_pupil_ridge,rsquared_grabs_face_ridge,rsquared_grabs_wheel_ridge);

% pad nan parcels up to 56 parcels
nan_parcels = nan(size(outputMat_grabs_regular,1),1);
outputMat_grabs_regular = [outputMat_grabs_regular(:,1:20) nan_parcels nan_parcels nan_parcels nan_parcels nan_parcels nan_parcels...
    outputMat_grabs_regular(:,21:end) nan_parcels nan_parcels nan_parcels nan_parcels];

% pad nan parcels up to 56 parcels
%nan_parcels = nan(size(outputMat_grabs_ridge,1),1);
%outputMat_grabs_ridge = [outputMat_grabs_ridge(:,1:20) nan_parcels nan_parcels nan_parcels nan_parcels nan_parcels nan_parcels...
%   outputMat_grabs_ridge(:,21:end) nan_parcels nan_parcels nan_parcels nan_parcels];
outputMat_grabs_ridge=[];
toc
else
    outputMat_grabs_regular=[];outputMat_grabs_ridge=[];
    disp(strcat(animal,num2str(k),'NO DATA'))
end
disp(strcat(animal,num2str(k),'done'))
end