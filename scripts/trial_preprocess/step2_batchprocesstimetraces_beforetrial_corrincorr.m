%% Concatenates trial data according to labels (correct, incorrect; running, not running; correct/running; 
%incorrect/running, correct/not running; incorrect/not running; correct/high pupil; incorrect/high pupil, correct/low pupil; incorrect/low pupil)
%over psychometric testing days and reshapes data into parcels*time.
clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
animals={'xt','xs','xx','xu','xz','xw'};
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(animal);
    [animalDate,dys,~,~]=animaltodays(animal);
    days_to_process=dys;
    CumulativeZscored=[];
    for dayy=1:length(unique(days_to_process)) %iterate over psychometric days
        res = load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\spt'), strcat(animal,num2str(days_to_process(dayy)),'imaging_time_traces_global.mat')));
        st=findClosestDouble(res.imaging_time_traces.t,-3);ed=findClosestDouble(res.imaging_time_traces.t,-0.1);
        CumulativeZscored = cat(3, CumulativeZscored, res.imaging_time_traces.Allen(:,st:ed,:)); %concatenate normalized trial data       
        clearvars res
    end
    psy2=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\', animalDate,'allen_dataproc_psychometric.mat')));
    Correct=reshape(CumulativeZscored(:,:,psy2.labels_trials.'==1),56,[]);
    Incorrect=reshape(CumulativeZscored(:,:,psy2.labels_trials.'==2),56,[]);
    %all running trials
    RunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.labels_trials.'<3),56,[]);
    %RunTrials=CumulativeTrialLabels(psy2.allRunningTrials==1&psy2.labels_trials.'<3,:);
    %all not running trials
    NotRunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.labels_trials.'<3),56,[]);
    %NotRunTrials=CumulativeTrialLabels(psy2.allRunningTrials==0&psy2.labels_trials.'<3,:);
    %high pupil trials
    HighPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'<3),56,[]);
    %HighPupTrials=CumulativeTrialLabels(psy2.allPupilTrials==2&psy2.labels_trials.'<3,:);
    %low pupil trials
    LowPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'<3),56,[]);
    %LowPupTrials=CumulativeTrialLabels(psy2.allPupilTrials==1&psy2.labels_trials.'<3,:);
    %all running, correct
    Correct_RunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.labels_trials.'==1),56,[]);
    %all running, incorrect
    Incorrect_RunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.labels_trials.'==2),56,[]);
    %all not running, correct
    Correct_NotRunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.labels_trials.'==1),56,[]);
    %all not running, incorrect
    Incorrect_NotRunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.labels_trials.'==2),56,[]);
    %high pupil, correct
    Correct_HighPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'==1),56,[]);
    %high pupil, incorrect
    Incorrect_HighPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'==2),56,[]);
    %low pupil, correct
    Correct_LowPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'==1),56,[]);
    %low pupil, incorrect
    Incorrect_LowPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'==2),56,[]);

    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'state_trials_networkanalysis'),'RunZscored', 'NotRunZscored','HighPupZscored','LowPupZscored','Correct_RunZscored','Incorrect_RunZscored','Correct_NotRunZscored','Incorrect_NotRunZscored','Correct','Incorrect','Correct_HighPupZscored','Incorrect_HighPupZscored','Correct_LowPupZscored','Incorrect_LowPupZscored');
    clearvars -except ir animals    
end