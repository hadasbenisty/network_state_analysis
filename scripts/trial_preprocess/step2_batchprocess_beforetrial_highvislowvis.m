%% for each animal, concatenate imaging data over days, normalize traces, and add trial labels
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
    HighVis=reshape(CumulativeZscored(:,:,psy2.allVisTrials==2&psy2.labels_trials.'<3),56,[]);
    LowVis=reshape(CumulativeZscored(:,:,psy2.allVisTrials==1&psy2.labels_trials.'<3),56,[]);
    %all running, HighVis
    HighVis_RunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.allVisTrials==2),56,[]);
    %all running, LowVis
    LowVis_RunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.allVisTrials==1),56,[]);
    %all not running, HighVis
    HighVis_NotRunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.allVisTrials==2),56,[]);
    %all not running, LowVis
    LowVis_NotRunZscored=reshape(CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.allVisTrials==1),56,[]);
    %high pupil, HighVis
    HighVis_HighPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.allVisTrials==2),56,[]);
    %high pupil, LowVis
    LowVis_HighPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.allVisTrials==1),56,[]);
    %low pupil, HighVis
    HighVis_LowPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.allVisTrials==2),56,[]);
    %low pupil, LowVis
    LowVis_LowPupZscored=reshape(CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.allVisTrials==1),56,[]);
    figure;hist(psy2.allVisTrials==2&psy2.labels_trials.'<3);title('high vis'); % to see how many high response trials were getting
    figure;hist(psy2.allVisTrials==1&psy2.labels_trials.'<3);title('lowvis');
    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'before_trials_networkanalysis'),'HighVis','LowVis','HighVis_RunZscored','LowVis_RunZscored','HighVis_NotRunZscored','LowVis_NotRunZscored','HighVis_HighPupZscored','LowVis_HighPupZscored','HighVis_LowPupZscored','LowVis_LowPupZscored');
    clearvars -except ir animals
end