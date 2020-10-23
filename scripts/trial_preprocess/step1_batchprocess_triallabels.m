%%code to (1) concatenate and normalize trial data over days of interest,
%(2) calculate slopes over a moving time window per trial, and
%(3)  create trial labels such as running, not running, high pupil,low
%pupil etc
clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'))
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions\libsvm-3.22\libsvm-3.22'))
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions\libsvm-3.22'))
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\')

%initialize animal and days to process
animal='xz';
[~,days_to_process,~,~] = animaltodays(animal);
mkdir(animal);
disp(animal)
parcel_method_i=2;labels_trials=[];labels_contrast=[];allWheels=[];allSlopes=[];allpupil=[];all_pupil_highthres=[];all_pupil_lowthres=[];all_slopes_highthres=[];all_slopes_lowthres=[];
%moving time window parameters
winStart = -2.5:0.1:3;
winEnd = winStart + 0.22;
%iterate over days, concatenating trials and slope values over trials
for dayy=1:length(days_to_process)
    animalDate=strcat(animal,num2str(days_to_process(dayy))); 
    animalTimeTrace=strcat(animalDate,'imaging_time_traces_global.mat');
    spike2 = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalDate,'spike2.mat')));
    wtcurrent=squeeze(spike2.spike2time_traces.wheelspeed);
    res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animalDate,'imaging_time_traces_global.mat')));
    pupil=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\videoFeatures'), strcat(animalDate,'pupil_trials.mat')));
    pupil_thres=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalDate,'running_ITI.mat')));
    alldataN = normByPreActivity(res.imaging_time_traces.t, res.imaging_time_traces.Allen, -3, -1); %zscore trial imaging data
    psy2=(fullfile(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\', animaltodays(animal),'allen_dataproc_psychometric.mat')));
    if exist(psy2,'file')==2 %if slope data already exists, load it instead of reprocessing
        loadallSlopes=1; 
        sl_info=alldataN;
    else
    loadallSlopes=0;
    sl_info=[];
    for t = 1:length(winStart) %iterate over time windows and store the slope by linear fit
        timeinds = find(res.imaging_time_traces.t>=winStart(t) &res.imaging_time_traces.t<=winEnd(t)); %time window
        for T=1:size(alldataN,3)
            for p = 1:size(alldataN,1)
                P = polyfit(timeinds,alldataN(p,timeinds, T),1); %slope calculation
                sl_info(t, p,T) = P(1);
            end
        end
    end
    allSlopes=cat(3, allSlopes, sl_info); 
    all_slopes_highthres=cat(2,all_slopes_highthres,ones(1,size(sl_info,3))*quantile(allSlopes(27,2,:),0.75)); %L-v1 slope quantile
    all_slopes_lowthres=cat(2,all_slopes_lowthres,ones(1,size(sl_info,3))*quantile(allSlopes(27,2,:),0.25));
    end
    allWheels=cat(2,allWheels,wtcurrent(:,1:(size(sl_info,3)))); %running speed 
    allpupil=cat(2,allpupil,pupil.pupil_time_trace.pupildata(:,1:(size(sl_info,3))));%pupil diameter
    all_pupil_highthres=cat(2,all_pupil_highthres,ones(1,size(sl_info,3))*pupil_thres.running_time_traces.zthres_High); %Pupil quantile
    all_pupil_lowthres=cat(2,all_pupil_lowthres,ones(1,size(sl_info,3))*pupil_thres.running_time_traces.zthres_Low);

    labels_trials = cat(1, labels_trials, res.trialslabels.blinksummary); 
    labels_contrast=cat(1, labels_contrast,res.trialslabels.contrastLabels);
    clearvars pupil alldataN res wtcurrent spike2 animalTimeTraceanimalDate pupil_thres
end

if loadallSlopes==1
    allSlopes=[];
    psy2info=load(psy2);
    allSlopes=psy2info.allSlopes;
    all_slopes_highthres=psy2info.all_slopes_highthres;
    all_slopes_lowthres=psy2info.all_slopes_lowthres;
    %plot high vs low visual slopes
    figure;
    hold on
    bar(cat(2,squeeze(mean(allSlopes(27,2,psy2info.allVisTrials==2&psy2info.labels_trials.'<3),3)),squeeze(mean(allSlopes(27,2,psy2info.allVisTrials==1&psy2info.labels_trials.'<3),3))));
    xlabel('High vs Low Vis');ylabel('Slopes');
    mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\',animalDate,'v1_slope_separation_corr_incorr'),'all');
else
end
%create running labels per trial

%divide trial into time window partitions
spike2 = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spike2Features'), strcat(animalDate,'spike2.mat')));
pupil=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\videoFeatures'), strcat(animalDate,'pupil_trials.mat')));

%iterate through trials, classifying those that have >= 1cm/s (putative) for
%any time window as running (1) , else as not running (0)
% e=findClosestDouble(pupil.pupil_time_trace.t,-6);
str=findClosestDouble(spike2.spike2time_traces.t,-6);
endr=findClosestDouble(spike2.spike2time_traces.t,1);

allRunningTrials=[];
for j=1:size(allWheels,2)
    Current=allWheels(str:endr,j);    
    abovethres=find(Current>10);
    belowthres=find(Current<10);
    if length(abovethres)==length(Current)
        trialinfo=1;
    elseif length(belowthres)==length(Current)
        trialinfo=0;
    else
        trialinfo=nan;
    end
    allRunningTrials(j)=trialinfo;
    clearvars abovethres belowthres trialinfo
end

%clearvars e m ml ml2 l

% e=findClosestDouble(pupil.pupil_time_trace.t,-6);
% m=findClosestDouble(pupil.pupil_time_trace.t,-3);
% ml=findClosestDouble(pupil.pupil_time_trace.t,0);
% l=findClosestDouble(pupil.pupil_time_trace.t,0.33);

%iterate through trials, classifying those that have >=0.6 quantile of
%pupil data as "high pupil",those with <=0.4 quantile of pupil
%data as low pupil, throw out everything between 0.4 and 0.6 quantile.
st=findClosestDouble(pupil.pupil_time_trace.t,-6);
ed=findClosestDouble(pupil.pupil_time_trace.t,0.33);

allPupilTrials=[];
for j=1:size(allWheels,2)
    trialinfo=[];
if allRunningTrials(j)==0
    Current=allpupil(st:ed,j);
    abovethres=find(Current>all_pupil_highthres(j));
    belowthres=find(Current<all_pupil_lowthres(j));
    if length(abovethres)==length(Current)
        trialinfo=2;
    elseif length(belowthres)==length(Current)
        trialinfo=1;
    else
        trialinfo=0;
    end
else
    trialinfo=0; 
end
    allPupilTrials(j)=trialinfo;
    clearvars Current trialinfo abovethres belowthres
end
disp(strcat(num2str(length(find(allPupilTrials==2))),'highpupil trials in', animal))
disp(strcat(num2str(length(find(allPupilTrials==1))),'lowpupil trials in', animal))
disp(strcat(num2str(length(find(allRunningTrials==1))),'running trials in', animal))
disp(strcat(num2str(length(find(allRunningTrials==0))),'not running trials in', animal))

%classify as high or low visual response trial from 100-320ms time window
%slope of l-v1
allVisTrials=[];
for j=1:size(allSlopes,3)
    trialinfo=[];
    if squeeze(allSlopes(27, 2, j))>=all_slopes_highthres(j)
        trialinfo=2;
    elseif squeeze(allSlopes(27, 2, j))<=all_slopes_lowthres(j)
        trialinfo=1;
    else
        trialinfo=0;
    end
    allVisTrials(j)=trialinfo;
    clearvars trialinfo
end

%sanity check plots
st=findClosestDouble(pupil.pupil_time_trace.t,-1.5);
ed=findClosestDouble(pupil.pupil_time_trace.t,1.5);
indx=st:ed;
close all;
figure;
hold on
shadedErrorBar(pupil.pupil_time_trace.t(indx),mean(allpupil(indx,  allPupilTrials==2),2),std(allpupil(indx,  allPupilTrials==2),0,2)/sqrt((size(allpupil(:, allPupilTrials==2),2)-1)),'lineprops','g');
hold on
shadedErrorBar(pupil.pupil_time_trace.t(indx),mean(allpupil(indx,  allPupilTrials==1),2),std(allpupil(indx,  allPupilTrials==1),0,2)/sqrt((size(allpupil(:, allPupilTrials==1),2)-1)),'lineprops','r');
xlabel('Time [Sec]');ylabel('Pupil Area');title('Pupil in High vs Low Pupil Trials')
mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\',animalDate,'v2_pupil_separation'),'all');

figure;
hold on
bar(cat(2,squeeze(mean(allSlopes(27,2,allVisTrials==2),3)),squeeze(mean(allSlopes(27,2,allVisTrials==1),3))));
xlabel('High vs Low Vis');ylabel('Slopes');
mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\',animalDate,'v2_v1_slope_separation'),'all');

%save workspace
save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\',animalDate,'v2_allen_dataproc_psychometric_workspace'));

%save select variabes
save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\',animalDate,'v2_allen_dataproc_psychometric'),...
    'allSlopes','allWheels','labels_trials','labels_contrast','allRunningTrials','winStart','winEnd','days_to_process','allPupilTrials','allpupil','allVisTrials','all_slopes_highthres','all_slopes_lowthres');

