%% Compute whether pupil diameter differs for correct vs incorrect low pupil trials
clear;
addpath(genpath('X:\Hadas\Meso-imaging\lan\results\code\Functions'));
cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude');
animals={'xt','xs','xx','xu','xz','xw'};
corr_pupil=[];incorr_pupil=[];corr_pupil_dist=[];incorr_pupil_dist=[];
for ir=1:length(animals)
    animal=char(animals(ir));
    mkdir(animal);
    [animalDate,dys,~,~]=animaltodays(animal);
    psy2=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\',animal,'\', animalDate,'allen_dataproc_psychometric.mat')));
    %correct vs incorret low pupil timestamps
    Correct_LowPup_Index=find(psy2.allPupilTrials==1&psy2.labels_trials.'==1);
    Incorrect_LowPup_Index=find(psy2.allPupilTrials==1&psy2.labels_trials.'==2);
    %get timestamps for -1.5s to 0.3s
    pupil=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\',animal,'psych\videoFeatures'), strcat(animalDate,'pupil_trials.mat')));
    st=findClosestDouble(pupil.pupil_time_trace.t,-1.5);
    ed=findClosestDouble(pupil.pupil_time_trace.t,0.3);    
    indx=st:ed;
    %concatenate pupil area from -1.5 to 0.3s in correct vs incorrect low pupil trials 
    %average over trials to get a timetrace
    corr_pupil=cat(2,corr_pupil,mean(psy2.allpupil(indx,Correct_LowPup_Index),2));
    incorr_pupil=cat(2,incorr_pupil,mean(psy2.allpupil(indx,Incorrect_LowPup_Index),2));
    %average over time, then average over trials to get average pupil area
    %for each animal per corr and incorr in high/low pupil
    corr_pupil_dist=cat(1,corr_pupil_dist,mean(mean(psy2.allpupil(indx,Correct_LowPup_Index),1)));
    incorr_pupil_dist=cat(1,incorr_pupil_dist,mean(mean(psy2.allpupil(indx,Incorrect_LowPup_Index),1)));
    clearvars -except ir animals corr_pupil incorr_pupil corr_pupil_dist incorr_pupil_dist
end
%paired, non parametric
p = signrank(corr_pupil_dist,incorr_pupil_dist)
figure;
%set(gcf,'renderer','Painters')
hold on
histogram(corr_pupil_dist,10,'facecolor','g','facealpha',.5,'edgecolor','none') %blue is incorrect
histogram(incorr_pupil_dist,10,'facecolor','r','facealpha',.5,'edgecolor','none')
xlabel('Average Pupil Area (t=-1.5s to 0.3s)')
title('Pupil Area in Correct (Green) vs Incorrect (Red)')
box off
axis tight
mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\basic_figs','pupildist','lowpupil_corr_incorr_dist'), 'all');


pupil=load(fullfile(strcat('X:\Hadas\Meso-imaging\lan\','xs','psych\videoFeatures'), strcat(animaltodays('xs'),'pupil_trials.mat')));
st=findClosestDouble(pupil.pupil_time_trace.t,-1.5);
ed=findClosestDouble(pupil.pupil_time_trace.t,0.3);
indx=st:ed;
figure;
%set(gcf,'renderer','Painters')
hold on
shadedErrorBar(pupil.pupil_time_trace.t(indx),mean(corr_pupil,2),std(corr_pupil,0,2)/sqrt((size(corr_pupil,2)-1)),'lineprops','g');
hold on
shadedErrorBar(pupil.pupil_time_trace.t(indx),mean(incorr_pupil,2),std(incorr_pupil,0,2)/sqrt((size(incorr_pupil,2)-1)),'lineprops','r');
xlabel('Time [Sec]');ylabel('Pupil Area');title('Pupil Area in Correct vs Incorrect Pupil Trials')
xlim([-1.5 0.3])
ylim([0 850])
mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\basic_figs','pupildist','lowpupil_corr_incorr_timetrace'), 'all');
