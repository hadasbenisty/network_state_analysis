function pre_process_trials_animals
%%Code to process the ITI period for all animals with extraction 
%of 3 arousal states: low pupil + quiescence, high pupil + quiescence, high pupil + locomotion


addpath(genpath('../../pre_processing'));
addpath(genpath('../../meta_data_processing'));
addpath(genpath('../../utils'));
animals={'xs','xx','xz','xw','xt','xu'};

for ai = 1:length(animals)
    extract_trials_imaging_by_state_loose(animals{ai})
%     extract_trials_imaging_by_state_super_strict(animals{ai});
end
isloose = true;
concatenateTrialsPeriodsByState(animals, isloose);

plot_time_spent(animals, isloose)

end
function extract_trials_imaging_by_state_loose(animalName)

params.fsimaging=10;%imaging sampling rate
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeBefore=4;% window of state before stim
params.Duration=4.4;% window length of state 
params.runningThSpeed = 1;% cm per sec
params.runningPercentTh = 10;% %of time animal is running on a trial so it would be considered a running trial
params.pupPercentTh= 10;% %of time animal is high pupil on a trial so it would be considered a high pupil trial
spike2pth = fullfile('X:\Hadas\Meso-imaging\lan\spike2data', animalName);
[~,days2process]=animaltodays(animalName);
disp(animalName);
%% for each mouse, load spont and airpuff folders and perform correlations
for day_i=1:length(days2process)
    disp(num2str(days2process(day_i)));
     load(fullfile(spike2pth, ['spike2data',animalName num2str(days2process(day_i)) '.mat']),'channels_data',...
        'timing', 'channels_data');
    
   
    [pupil_time, pupil_Norm] = load_pupil_data(animalName, days2process(day_i), channels_data.pupil_frame, params.fsspike2);
    wheel_speed = channels_data.wheelspeed;
    wheel_time = (1:length(wheel_speed))/params.fsspike2;
    stimtimes = timing.stimstart/params.fsspike2;
    stimtimes=stimtimes(1:75);
    before_time = stimtimes - params.TimeBefore;
    for stim_i = 1:length(stimtimes)
        ind1 = findClosestDouble(before_time(stim_i), wheel_time);
        ind2 = ind1 + params.Duration*params.fsspike2;
        wheel_trials(:,stim_i) = wheel_speed(ind1:ind2);        
        
        ind1 = findClosestDouble(before_time(stim_i), pupil_time);
        ind2 = ind1 + params.Duration*params.fspupilcam;
        pupil_trials(:,stim_i) = pupil_Norm(ind1:ind2);
    end
    trial_label = nan(length(stimtimes),1);
    for stim_i = 1:length(stimtimes)
        running_times = wheel_trials(:,stim_i) > params.runningThSpeed;
        if 100*sum(running_times)/length(wheel_trials(:,stim_i)) > params.runningPercentTh
           trial_label(stim_i) = 3;
        end
    end
    q_pupil_vals = pupil_trials(:, trial_label~=1);
    zthres_High=quantile(q_pupil_vals(:),0.60);
    for stim_i = 1:length(stimtimes)
        if isnan(trial_label(stim_i))
            highpup_times = pupil_trials(:,stim_i) > zthres_High;
            if 100*sum(highpup_times)/length(pupil_trials(:,stim_i)) > params.pupPercentTh
                trial_label(stim_i) = 2;
            else
                trial_label(stim_i) = 1;
            end
        end
    end
    
  
    trials_labels_arousal_pup_loose = trial_label;
    
    trials_labels_arousal_pup_loose_lut = {'pupil_low_q', 'pupil_high_q','pupil_high_loc'};
    if any(isnan(trials_labels_arousal_pup_loose))
        error('there is at least one trial without a label')
    end
    
    
    
    
    
       datapath = ['X:\Hadas\Meso-imaging\lan\' animalName 'psych\spt\'];

   
          load(fullfile(datapath,[animalName num2str(days2process(day_i)) 'imaging_time_traces_global.mat']),...
              'trials_labels_arousal_pup_lut','trials_labels_arousal_facemap_lut','trialslabels','roiLabelsbyAllen','maskByAllen','regionLabel','isLeftLabel','imaging_time_traces');
    trialslabels.trials_labels_arousal_pup_loose = trials_labels_arousal_pup_loose;
    
    save(fullfile(datapath,[animalName num2str(days2process(day_i)) 'imaging_time_traces_global.mat']),...
        'trialslabels','roiLabelsbyAllen','maskByAllen','regionLabel','isLeftLabel','imaging_time_traces',...
        'trials_labels_arousal_pup_lut','trials_labels_arousal_facemap_lut', 'trials_labels_arousal_pup_loose_lut')
end

end

function plot_time_spent(animals, isloose)

if isloose
    loosestr = 'loose';
else
    loosestr = '';
end

for ai = 1:length(animals)
    load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animals{ai},'\',animals{ai},'trials_3states', loosestr),...
        'low_pup_q','high_pup_q','high_pup_l');
    Nall(ai, 1) = sum(low_pup_q.trialslabels.blinksummary<3);
    Ncor(ai, 1) = sum(low_pup_q.trialslabels.blinksummary==1);
    Ninc(ai, 1) = sum(low_pup_q.trialslabels.blinksummary==2);
    
    Nall(ai, 2) = sum(high_pup_q.trialslabels.blinksummary<3);
    Ncor(ai, 2) = sum(high_pup_q.trialslabels.blinksummary==1);
    Ninc(ai, 2) = sum(high_pup_q.trialslabels.blinksummary==2);
    
    Nall(ai, 3) = sum(high_pup_l.trialslabels.blinksummary<3);
    Ncor(ai, 3) = sum(high_pup_l.trialslabels.blinksummary==1);
    Ninc(ai, 3) = sum(high_pup_l.trialslabels.blinksummary==2);
end
n=length(animals);
Mall = nanmean(Nall);    
Sall = nanstd(Nall)/sqrt(n-1);  
Mcor = nanmean(Ncor);  
Scor = nanstd(Ncor)/sqrt(n-1);  
Minc = nanmean(Ninc);  
Sinc = nanstd(Ninc)/sqrt(n-1);  

figure;subplot(3,1,1);
barwitherr(Sall, Mall);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('All Trials');
subplot(3,1,2);
barwitherr(Scor, Mcor);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('Correct Trials');
subplot(3,1,3);
barwitherr(Sinc, Minc);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('Incorrect Trials');
mysave(gcf,['X:\Lav\ProcessingDirectory\parcor_undirected\time_spent_bytrials_absolute_numbers' loosestr]) 



n=length(animals);
Mall = nanmean(Nall./sum(Nall,2));    
Sall = nanstd(Nall./sum(Nall,2))/sqrt(n-1);  
Mcor = nanmean(Ncor./sum(Ncor,2));  
Scor = nanstd(Ncor./sum(Ncor,2))/sqrt(n-1);  
Minc = nanmean(Ninc./sum(Ninc,2));  
Sinc = nanstd(Ninc./sum(Ninc,2))/sqrt(n-1);  

figure;subplot(4,1,1);
barwitherr(Sall, Mall);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('All Trials');
subplot(4,1,2);
barwitherr(Scor, Mcor);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('Correct Trials');
subplot(4,1,3);
barwitherr(Sinc, Minc);
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('Incorrect Trials');
%% 
n=length(animals);
Nall = Ncor + Ninc;
Mcor = nanmean(Ncor./Nall);  
Scor = nanstd(Ncor./Nall)/sqrt(n-1);  
Minc = nanmean(Ninc./Nall);  
Sinc = nanstd(Ninc./Nall)/sqrt(n-1);  

subplot(4,1,4);

barwitherr([Scor; Sinc]' , [Mcor;Minc]');
set(gca,'XTickLabel', {'low pup q','high pup q','high pup l'});
title('Time Spent By State');legend('Correct', 'Incorrect');


mysave(gcf,['X:\Lav\ProcessingDirectory\parcor_undirected\time_spent_bytrials' loosestr]) 



end
function concatenateTrialsPeriodsByState(animals, isloose)
%% Concatenates spontaneous state “trials” from the previous step over all days in psyc testing,
%reshaped to be parcels over time (for running not running).
fltstr='spt';
addpath(genpath('../meta_data_processing'))
respath = 'X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory';
for ir=1:length(animals)
    animal=char(animals(ir));
    datapath = ['X:\Hadas\Meso-imaging\lan\' animal 'psych\' fltstr '\'];
    disp(animal)
    
    
    mkNewDir(fullfile(respath, 'allen_Slope_Amplitude',animal));
    [~,days_to_process]=animaltodays(animal);
    days_to_process = unique(days_to_process);
    
    low_pup_q.imaging_time_traces=[];
    low_pup_q.trialslabels.blinksummary=[];
    low_pup_q.trialslabels.contrastLabels = [];
    
    high_pup_q.imaging_time_traces=[];
    high_pup_q.trialslabels.blinksummary=[];
    high_pup_q.trialslabels.contrastLabels = [];
    
    high_pup_l.imaging_time_traces=[];
    high_pup_l.trialslabels.blinksummary=[];
    high_pup_l.trialslabels.contrastLabels = [];
    
    if isloose
        lutvar = 'trials_labels_arousal_pup_loose_lut';
        labelsvar = 'trials_labels_arousal_pup_loose';
    else
        lutvar = 'trials_labels_arousal_pup_lut';
        labelsvar = 'trials_labels_arousal_pup';
    end
    for dayy=1:length(days_to_process) %iterate over psychometric days
        disp(days_to_process(dayy))
         datafile = fullfile(datapath,[animal num2str(days_to_process(dayy)) 'imaging_time_traces_global.mat']);
    
        if exist(datafile,'file')
             load(datafile, 'trialslabels','imaging_time_traces',...
        lutvar);
 
    
            % state 1: low pup + q   
            low_pup_q = get_trials_by_state(imaging_time_traces, trialslabels, labelsvar, 'pupil_low_q', eval(lutvar), low_pup_q);
            % state 2: high pup + q  
            high_pup_q = get_trials_by_state(imaging_time_traces, trialslabels, labelsvar, 'pupil_high_q', eval(lutvar), high_pup_q);
            % state 3: high pup + loc       
            high_pup_l = get_trials_by_state(imaging_time_traces, trialslabels, labelsvar, 'pupil_high_loc', eval(lutvar), high_pup_l);
            
        end
    end
    t = imaging_time_traces.t;
     if isloose
         loosestr = 'loose';
     else
         loosestr = '';
     end
    save(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states', loosestr),...
        'low_pup_q','high_pup_q','high_pup_l','days_to_process', 't');
end


end
function extract_trials_imaging_by_state_super_strict(animalName)

params.fsimaging=10;%imaging sampling rate
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=5;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=5;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=3;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal)
params.minSitDuration=3;%minimum sit duration during quiescnece state

spike2pth = fullfile('X:\Hadas\Meso-imaging\lan\spike2data', animalName);
[~,days2process]=animaltodays(animalName);
disp(animalName);
%% for each mouse, load spont and airpuff folders and perform correlations
for day_i=1:length(days2process)
    disp(num2str(days2process(day_i)));
     load(fullfile(spike2pth, ['spike2data',animalName num2str(days2process(day_i)) '.mat']),'channels_data',...
        'timing', 't_imaging');
    %     X:\Hadas\Meso-imaging\lan\facemap\Current_Processing\xu\xu_D20_proc.mat
    
    facemap_data_file=dir(strcat('X:\Hadas\Meso-imaging\lan\facemap\Current_Processing\', animalName, ['\*' animalName '_D' ...
    num2str(days2process(day_i)) '*_proc.mat']));
    facemapfile=(fullfile(strcat('X:\Hadas\Meso-imaging\lan\facemap\Current_Processing\', animalName),facemap_data_file.name));    
    %facemapfile = ['X:\Hadas\Meso-imaging\lan\facemap\Current_Processing\' animalName ...
    %   '\' animalName '_D' num2str(days2process(day_i)) '_proc.mat'];
    if exist(facemapfile, 'file')
        q=load(facemapfile);
        face_Norm = q.proc.motSVD{1,1}(:,1);
    else
        face_Norm = [];
    end
    [pupil_time, pupil_Norm] = load_pupil_data(animalName, days2process(day_i), channels_data.pupil_frame, params.fsspike2);
    wheel_speed = channels_data.wheelspeed;
    wheel_time = (1:length(wheel_speed))/params.fsspike2;
    
    %% get locomotion on/off and quiescence on off times
    %locomotion periods should be at least some criterion s long with some criterion s since locomotion
    %onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)
    raw_spike2 = fullfile('X:\Lan\Meso-imaging\', animalName, [animalName '_D' ...
        num2str(days2process(day_i))  '.mat']);
    raw_channels = load(raw_spike2,  'data','timestamp'); %load raw wheel data from lan because we need the rotations for
    
    
    [~, ~, ~, wheelOn, wheelOff] = get_wheel_on_off(raw_channels.data(:,5).', params.fsspike2, raw_channels.timestamp, timing);
    
    %find wheel on and wheel off times during imaging period only
    minRunDur=params.minRunDuration+params.TimeSinceLocOn+params.TimeBeforeLocOff; %minimum actual locomotion duration including time since locomotion onset, time before locomotion offset and the minimum time period for data analysis
    idx=wheelOn<(t_imaging(end)) & wheelOff>(t_imaging(1));
    wheelOn_t1=wheelOn(idx);
    wheelOff_t1=wheelOff(idx);
    
    for whe=1:length(wheelOn_t1)
        if wheelOn_t1(whe)<t_imaging(end) && wheelOff_t1(whe)>t_imaging(end)
            wheelOff_t1(whe)=t_imaging(end)-1;%if locomotion starts before end of imaging but continues after, only extract state until imaging time end minus a second
        end
        
        if wheelOn_t1(whe)<t_imaging(1) && wheelOff_t1(whe)>t_imaging(1)
            wheelOn_t1(whe)=t_imaging(1)+1;%if locomotion starts before start of imaging but continues after, only extract state from imaging time start plus a second
        end
    end
    
    wheelOn_t1=(wheelOn_t1(:))';
    wheelOff_t1=(wheelOff_t1(:))';
    
    %find wheel on off  times when airpuffs are not given
    if isempty(wheelOn_t1)
        disp('no wheel on');
        continue;
    end
    allEvts=sort(timing.stimstart/params.fsspike2,'ascend');
    allEvts=allEvts(:)';
    allEvtsPre=allEvts-params.TimeSinceEvent;
    allEvtsPost=allEvts+params.TimeSinceEvent;
    
    newwheelOn=sort([wheelOn_t1,allEvtsPre,allEvtsPost],'ascend');
    newwheelOff=sort([wheelOff_t1,allEvtsPre,allEvtsPost],'ascend');
    wheelOn_int = wheelOn_t1;
    wheelOff_int = wheelOff_t1;
%     Index=nan(1,length(newwheelOn));
%     for r=1:length(newwheelOn)
%         tmp1=newwheelOn(r);
%         tmp2=newwheelOff(r);
%         
%         tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
%         if sum(tmp3)==0
%             Index(r)=r;
%         end
%     end
%     wheelOn_int=newwheelOn(Index(~isnan(Index)));
%     wheelOff_int=newwheelOff(Index(~isnan(Index)));
    
    %makes sure the state is at least as long as the minimum run duration
    idx1=find((wheelOff_int-wheelOn_int)>=(minRunDur));
    wheelOn_int1=wheelOn_int(idx1);
    wheelOff_int1=wheelOff_int(idx1);
    
    %finalize the times to get sustained state only
    wheelOn_final=wheelOn_int1+params.TimeSinceLocOn;
    wheelOff_final=wheelOff_int1-params.TimeBeforeLocOff;
    
    %% queiscence should be at least some criterion s long with some criterion s since locomotion offset
    %and some criterion s before subsequent locomotion onset,excluding any events (airpuff/stim)
    sitOn=[0;wheelOff(1:end-1)]; %use 0 as the first sit on time;
    sitOff=wheelOn;%use wheelOn times as sit off times;
    
    %find sit on and sit off times during imaging period only
    minSitDur=params.minSitDuration+params.TimeSinceSitOn+params.TimeBeforeSitOff; %actual minimum sit duration accouting for the onset time, offset time and minimum duration of the sustained quiescence epoch  used for analysis
    
    idx=sitOn<(t_imaging(end)) & sitOff>(t_imaging(1));
    sitOn_t1=sitOn(idx);
    sitOff_t1=sitOff(idx);
    
    for whe=1:length(sitOn_t1)
        if sitOn_t1(whe)<t_imaging(end) && sitOff_t1(whe)>t_imaging(end)
            sitOff_t1(whe)=t_imaging(end)-1;%if queiscence starts before end of imaging but continues after, only extract state until imaging time end minus a second
        end
        
        if sitOn_t1(whe)<t_imaging(1) && sitOff_t1(whe)>t_imaging(1)
            sitOn_t1(whe)=t_imaging(1)+1;%if queiscence starts before start of imaging but continues after, only extract state from imaging time start plus a second
        end
    end
    
    %remove any quiescence period where mouse's speed is above 0.03 m/s
    % because sometimes these epcohs,especially short ones, get missed by locomotion changepoint algorithm
    speedThres=0.03;
    tmpOn=sitOn_t1; tmpOff=sitOff_t1;
    for rr=1:length(sitOn_t1)
        currspeed = wheel_speed(wheel_time>sitOn_t1(rr)& wheel_time<sitOff_t1(rr));
        highSpeedIdx=abs(currspeed)>speedThres;
        if ~isempty(highSpeedIdx)
            firstIdx=wheel_time(find(wheel_time>sitOn_t1(rr),1)+highSpeedIdx(1))-0.1;
            lastIdx=wheel_time(find(wheel_time>sitOn_t1(rr),1)+highSpeedIdx(end))+0.1;
            tmpOff(rr)=firstIdx;
            tmpOff(end+1)=sitOff_t1(rr);
            tmpOn(end+1)=lastIdx;
        end
    end
    tmpOn=sort(tmpOn,'ascend'); tmpOff=sort(tmpOff,'ascend');
    
    sitOn_t1=(tmpOn(:))';
    sitOff_t1=(tmpOff(:))';
    sitOn_int=sitOn_t1;
    sitOff_int=sitOff_t1;
    
    %find sit on off  times when airpuffs are not given
%     allEvts=sort(timing.stimstart/params.fsspike2,'ascend');
%     allEvts=allEvts(:)';
%     allEvtsPre=allEvts-params.TimeSinceEvent;
%     allEvtsPost=allEvts+params.TimeSinceEvent;
%     
%     newSitOn=sort([sitOn_t1,allEvtsPre,allEvtsPost],'ascend');
%     newSitOff=sort([sitOff_t1,allEvtsPre,allEvtsPost],'ascend');
%     
%     Index=nan(1,length(newSitOn));
%     for r=1:length(newSitOn)
%         tmp1=newSitOn(r);
%         tmp2=newSitOff(r);
%         
%         tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
%         if sum(tmp3)==0
%             Index(r)=r;
%         end
%     end
%     sitOn_int=newSitOn(Index(~isnan(Index)));
%     sitOff_int=newSitOff(Index(~isnan(Index)));
    
%     
    %makes sure state is at least as long as the minimum sit duration
    idx1=find((sitOff_int-sitOn_int)>=(minSitDur));
    sitOn_int1=sitOn_int(idx1);
    sitOff_int1=sitOff_int(idx1);
    
    %finalize the times to get sustained state only
    sitOn_final=sitOn_int1+params.TimeSinceSitOn;
    sitOff_final=sitOff_int1-params.TimeBeforeSitOff;
    
    %% do change point detection on pupil and face to get pupil high/low arousal or face high/low movement times during sustained quiescence state
    %get Z-thresholds based on pupil data during quiescence, when mouse
    %isn't moving and when aripuffs are not given
    b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
    pupilTime_Idx=cell(1,length(sitOn_int));
    for st=1:length(sitOn_int)
        pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
    end
    pupilTime_quiescence=cell2mat(pupilTime_Idx');
    pupil_qui_times=pupil_time(pupilTime_quiescence);
    pupil_quiescence=pupil_Norm(pupilTime_quiescence);
    zthres_High=quantile(pupil_quiescence,0.60);
    zthres_Low=quantile(pupil_quiescence,0.40);
    
    %get on and off timestamps for high and low arousal based on pupil
    [h1,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp ] =changepoints1(pupil_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
    [h2,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints1(-pupil_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
    %get z thresholds based on face data during quiescence only
    if ~isempty(face_Norm)
        b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
        face_quiescence=face_Norm(pupilTime_quiescence);
        zthres_High=quantile(face_quiescence,0.60);
        zthres_Low=quantile(face_quiescence,0.40);
        
        %get on and off timestamps for high and low face movment
        [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints1(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
        [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints1(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
    end
    % determine that pupil and face high/low arousal/movement time are at least minimum criterion seconds long
    idx1=find((Pupil_HighArousal_OffTStamp-Pupil_HighArousal_OnTStamp)>=params.minArousalDuration);
    Pupil_HighArousal_OnT_int=Pupil_HighArousal_OnTStamp(idx1); Pupil_HighArousal_OffT_int=Pupil_HighArousal_OffTStamp(idx1);
    idx2=find((Pupil_LowArousal_OffTStamp-Pupil_LowArousal_OnTStamp)>=params.minArousalDuration);
    Pupil_LowArousal_OnT_int=Pupil_LowArousal_OnTStamp(idx2); Pupil_LowArousal_OffT_int=Pupil_LowArousal_OffTStamp(idx2);
    if ~isempty(face_Norm)
        if length(Face_LowArousal_OffTStamp)<length(Face_LowArousal_OnTStamp)
            Face_LowArousal_OnTStamp=Face_LowArousal_OnTStamp(1:length(Face_LowArousal_OffTStamp));
        end
        idx3=find((Face_HighArousal_OffTStamp-Face_HighArousal_OnTStamp)>=params.minArousalDuration);
        Face_HighArousal_OnT_int=Face_HighArousal_OnTStamp(idx3); Face_HighArousal_OffT_int=Face_HighArousal_OffTStamp(idx3);
        idx4=find((Face_LowArousal_OffTStamp-Face_LowArousal_OnTStamp)>=params.minArousalDuration);
        Face_LowArousal_OnT_int=Face_LowArousal_OnTStamp(idx4); Face_LowArousal_OffT_int=Face_LowArousal_OffTStamp(idx4);
    end
    %% pupil on/off for q
    % get pupil on and face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
    toDelete=ones(1,length(Pupil_HighArousal_OnT_int));
    for rj=1:length(Pupil_HighArousal_OnT_int)
        tmp = find (Pupil_HighArousal_OnT_int(rj)>=sitOn_final & Pupil_HighArousal_OffT_int(rj)<=sitOff_final);
        toDelete(rj)=isempty(tmp);
    end
    Pupil_HighArousal_On_final_qu=Pupil_HighArousal_OnT_int(~toDelete);
    Pupil_HighArousal_Off_final_qu=Pupil_HighArousal_OffT_int(~toDelete);
    
    toDelete=ones(1,length(Pupil_LowArousal_OnT_int));
    for rj=1:length(Pupil_LowArousal_OnT_int)
        tmp = find (Pupil_LowArousal_OnT_int(rj)>=sitOn_final & Pupil_LowArousal_OffT_int(rj)<=sitOff_final);
        toDelete(rj)=isempty(tmp);
    end
    Pupil_LowArousal_On_final_qu=Pupil_LowArousal_OnT_int(~toDelete);
    Pupil_LowArousal_Off_final_qu=Pupil_LowArousal_OffT_int(~toDelete);
    %% pupil on/off for locomotion
    % get pupil on and face on times if both on and off times occur entirely
    %during sustained locomotion states identified in the previous step
    [Pupil_HighArousal_On_final_loc, Pupil_HighArousal_Off_final_loc,...
        Pupil_LowArousal_On_final_loc, Pupil_LowArousal_Off_final_loc] = ...
        getHighLowonOf4Locomotion(wheelOn_final, wheelOff_final, t_imaging, ...
        Pupil_HighArousal_OnT_int, Pupil_HighArousal_OffT_int, ...
        Pupil_LowArousal_OnT_int, Pupil_LowArousal_OffT_int);
    
    if ~isempty(face_Norm)
    %% facemap on/off for locomotion
    % get pupil on and face on times if both on and off times occur entirely
    %during sustained locomotion states identified in the previous step
    [Face_HighArousal_On_final_loc,Face_HighArousal_Off_final_loc,...
        Face_LowArousal_On_final_loc, Face_LowArousal_Off_final_loc] = ...
        getHighLowonOf4Locomotion(wheelOn_final, wheelOff_final, t_imaging, ...
        Face_HighArousal_OnT_int, Face_HighArousal_OffT_int, ...
        Face_LowArousal_OnT_int, Face_LowArousal_OffT_int);
    
    
    
    %% facemap
    
        toDelete=ones(1,length(Face_HighArousal_OnT_int));
        for rj=1:length(Face_HighArousal_OnT_int)
            tmp = find (Face_HighArousal_OnT_int(rj)>=sitOn_final & Face_HighArousal_OffT_int(rj)<=sitOff_final);
            toDelete(rj)=isempty(tmp);
        end
        Face_HighArousal_On_final_q=Face_HighArousal_OnT_int(~toDelete);
        Face_HighArousal_Off_final_q=Face_HighArousal_OffT_int(~toDelete);
        
        toDelete=ones(1,length(Face_LowArousal_OnT_int));
        for rj=1:length(Face_LowArousal_OnT_int)
            tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
            toDelete(rj)=isempty(tmp);
        end
        Face_LowArousal_On_final_q=Face_LowArousal_OnT_int(~toDelete);
        Face_LowArousal_Off_final_q=Face_LowArousal_OffT_int(~toDelete);
    end
    %% plot on and off times for sanity check
    %plot locomotion state on off times
    h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
    set(0,'CurrentFigure',h);ax1=subplot(6,1,1);
    plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('LocomotionState');
    for tt=1:length(wheelOn_final)
        plot([wheelOn_final(tt),wheelOn_final(tt)], ylimits, 'g');
        plot([wheelOff_final(tt),wheelOff_final(tt)], ylimits, 'r');
    end
    %plot airpuff timea andn imaging onset/offset
    if ~isempty(allEvts)
        for tt=1:length(allEvts)
            plot([allEvts(tt),allEvts(tt)], ylimits, 'k');
        end
    end
    plot([t_imaging(1),t_imaging(1)], ylimits, 'm');
    plot([t_imaging(end),t_imaging(end)], ylimits, 'm');
    
    %plot quiescence state on off times
    set(0,'CurrentFigure',h);ax2=subplot(6,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
    for tt=1:length(sitOn_final)
        plot([sitOn_final(tt),sitOn_final(tt)], ylimits, 'g');
        plot([sitOff_final(tt),sitOff_final(tt)], ylimits, 'r');
    end
    %plot airpuff timea and imaging onset/offset times
    if exist('allEvts','var')
        for tt=1:length(allEvts)
            plot([allEvts(tt),allEvts(tt)], ylimits, 'k');
        end
    end
    plot([t_imaging(1),t_imaging(1)], ylimits, 'm');
    plot([t_imaging(end),t_imaging(end)], ylimits, 'm');
    
    %plot pupil and face state on off times
    set(0,'CurrentFigure',h);ax3=subplot(6,1,3);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilHighState+q');
    for tt=1:length(Pupil_HighArousal_On_final_qu)
        plot([Pupil_HighArousal_On_final_qu(tt),Pupil_HighArousal_On_final_qu(tt)], ylimits, 'g');
        plot([Pupil_HighArousal_Off_final_qu(tt),Pupil_HighArousal_Off_final_qu(tt)], ylimits, 'r');
    end
    
    set(0,'CurrentFigure',h);ax4=subplot(6,1,4);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilLowState+q');
    for tt=1:length(Pupil_LowArousal_On_final_qu)
        plot([Pupil_LowArousal_On_final_qu(tt),Pupil_LowArousal_On_final_qu(tt)], ylimits, 'g');
        plot([Pupil_LowArousal_Off_final_qu(tt),Pupil_LowArousal_Off_final_qu(tt)], ylimits, 'r');
    end
    %plot pupil and face state on off times
    set(0,'CurrentFigure',h);ax5=subplot(6,1,5);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilHighState+loc');
    for tt=1:length(Pupil_HighArousal_On_final_loc)
        plot([Pupil_HighArousal_On_final_loc(tt),Pupil_HighArousal_On_final_loc(tt)], ylimits, 'g');
        plot([Pupil_HighArousal_Off_final_loc(tt),Pupil_HighArousal_Off_final_loc(tt)], ylimits, 'r');
    end
    
    set(0,'CurrentFigure',h);ax6=subplot(6,1,6);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilLowState+loc');
    for tt=1:length(Pupil_LowArousal_On_final_loc)
        plot([Pupil_LowArousal_On_final_loc(tt),Pupil_LowArousal_On_final_loc(tt)], ylimits, 'g');
        plot([Pupil_LowArousal_Off_final_loc(tt),Pupil_LowArousal_Off_final_loc(tt)], ylimits, 'r');
    end
    
    if ~isempty(face_Norm) && 0
        set(0,'CurrentFigure',h);ax5=subplot(6,1,5);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceHighState');
        for tt=1:length(Face_HighArousal_On_final)
            plot([Face_HighArousal_On_final(tt),Face_HighArousal_On_final(tt)], ylimits, 'g');
            plot([Face_HighArousal_Off_final(tt),Face_HighArousal_Off_final(tt)], ylimits, 'r');
        end
        
        set(0,'CurrentFigure',h);ax6=subplot(6,1,6);plot(pupil_time,face_Norm);ylimits=ylim; hold on; title('FaceLowState');
        for tt=1:length(Face_LowArousal_On_final)
            plot([Face_LowArousal_On_final(tt),Face_LowArousal_On_final(tt)], ylimits, 'g');
            plot([Face_LowArousal_Off_final(tt),Face_LowArousal_Off_final(tt)], ylimits, 'r');
        end
        
    end
    linkaxes([ax1, ax2, ax3,ax4,ax5,ax6],'x');
    
    
    allEvts=sort(timing.stimstart/params.fsspike2,'ascend');

    for ti = 1:length(allEvts)
%         is_running_trial(ti) = any(allEvts(ti) - wheelOn_final > 0 & allEvts(ti) - wheelOff_final < 0);
        is_puplow_on_q(ti) = any(allEvts(ti) - Pupil_LowArousal_On_final_qu > 3 & allEvts(ti) - Pupil_LowArousal_Off_final_qu < -1);
        is_puphigh_on_q(ti) = any(allEvts(ti) - Pupil_HighArousal_On_final_qu > 3 & allEvts(ti) - Pupil_HighArousal_Off_final_qu < -1);
        is_puplow_on_loc(ti) = any(allEvts(ti) - Pupil_LowArousal_On_final_loc > 3 & allEvts(ti) - Pupil_LowArousal_Off_final_loc < -1);
        is_puphigh_on_loc(ti) = any(allEvts(ti) - Pupil_HighArousal_On_final_loc > 3 & allEvts(ti) - Pupil_HighArousal_Off_final_loc < -1);
        is_low_face_on_q(ti) =  any(allEvts(ti) - Face_LowArousal_On_final_q > 3 & allEvts(ti) - Face_LowArousal_Off_final_q < -1);
        is_high_face_on_q(ti) =  any(allEvts(ti) - Face_HighArousal_On_final_q > 3 & allEvts(ti) - Face_HighArousal_On_final_q < -1);
       is_low_face_on_loc(ti) =  any(allEvts(ti) - Face_LowArousal_On_final_loc > 3 & allEvts(ti) - Face_LowArousal_Off_final_loc < -1);
       is_high_face_on_loc(ti) =  any(allEvts(ti) - Face_HighArousal_On_final_loc > 3 & allEvts(ti) - Face_HighArousal_Off_final_loc < -1);
        
    end
    if any(is_puplow_on_q + is_puphigh_on_q + is_puplow_on_loc + is_puphigh_on_loc>1)
        error('trials cannot have more than one state')
    end
    trials_labels_arousal_pup = int16(is_puplow_on_q);
    trials_labels_arousal_pup(is_puphigh_on_q == 1) = 2;
    trials_labels_arousal_pup(is_puplow_on_loc == 1) = 3;
    trials_labels_arousal_pup(is_puphigh_on_loc == 1) = 4;
    trials_labels_arousal_pup_lut = {'pupil_low_q', 'pupil_high_q','pupil_low_loc','pupil_high_loc'};
    if any(is_low_face_on_q + is_high_face_on_q + is_low_face_on_loc + is_high_face_on_loc>1)
        error('trials cannot have more than one state')
    end
    trials_labels_arousal_facemap = int16(is_low_face_on_q);
    trials_labels_arousal_facemap(is_high_face_on_q == 1) = 2;
    trials_labels_arousal_facemap(is_low_face_on_loc == 1) = 3;
    trials_labels_arousal_facemap(is_high_face_on_loc == 1) = 4;
    
        trials_labels_arousal_facemap_lut = {'face_low_q', 'face_high_q','face_low_loc','face_high_loc'};
 
    
    close all;
    
    
       datapath = ['X:\Hadas\Meso-imaging\lan\' animalName 'psych\spt\'];

   
          load(fullfile(datapath,[animalName num2str(days2process(day_i)) 'imaging_time_traces_global.mat']),...
              'trialslabels','roiLabelsbyAllen','maskByAllen','regionLabel','isLeftLabel','imaging_time_traces');
    trialslabels.trials_labels_arousal_pup = trials_labels_arousal_pup;
    trialslabels.trials_labels_arousal_facemap = trials_labels_arousal_facemap;
    
    save(fullfile(datapath,[animalName num2str(days2process(day_i)) 'imaging_time_traces_global.mat']),...
        'trialslabels','roiLabelsbyAllen','maskByAllen','regionLabel','isLeftLabel','imaging_time_traces',...
        'trials_labels_arousal_pup_lut','trials_labels_arousal_facemap_lut')
end

end

function trials_data = get_trials_by_state(imaging_time_traces, trialslabels, labelsname, statestr, trials_labels_arousal_pup_lut, trials_data)

state_label = find(strcmp(trials_labels_arousal_pup_lut, statestr));
inds = trialslabels.(labelsname) == state_label;
inds=inds(1:size(imaging_time_traces.Allen,3));
trials_data.imaging_time_traces=cat(3,trials_data.imaging_time_traces,imaging_time_traces.Allen(:,:,inds));
trials_data.trialslabels.blinksummary = cat(1, trials_data.trialslabels.blinksummary, trialslabels.blinksummary(inds));
trials_data.trialslabels.contrastLabels = cat(1, trials_data.trialslabels.contrastLabels, trialslabels.contrastLabels(inds));
end