clear;
addpath(genpath('X:\Lav\network_state_analysis\functions'));
animals={'xs','xt','xu','xx','xz','xw'};
cd('X:\Lav\network_state_analysis\scripts\permute_parcor_undirected\spon_highpup_lowpup');
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
    %     x.Correct=CumulativeZscored(:,:,psy2.labels_trials.'==1);
    %     x.Incorrect=CumulativeZscored(:,:,psy2.labels_trials.'==2);
    %     x.HighPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'<3);
    %     x.LowPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'<3);
    x.Correct_RunZscored=CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.labels_trials.'==1);
    x.Incorrect_RunZscored=CumulativeZscored(:,:,psy2.allRunningTrials==1&psy2.labels_trials.'==2);
    x.Correct_NotRunZscored=CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.labels_trials.'==1);
    x.Incorrect_NotRunZscored=CumulativeZscored(:,:,psy2.allRunningTrials==0&psy2.labels_trials.'==2);
    
    x.Correct_HighPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'==1);
    x.Incorrect_HighPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==2&psy2.labels_trials.'==2);
    x.Correct_LowPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'==1);
    x.Incorrect_LowPupZscored=CumulativeZscored(:,:,psy2.allPupilTrials==1&psy2.labels_trials.'==2);
    for ki=1:100
        permute_corrincorr_states(x,ki,animal);
    end
    clearvars -except ir animals
end