function pre_process_spont_animals_crispr
%%Code to process the ITI period for all animals with extraction
%of 3 arousal states: low pupil + quiescence, high pupil + quiescence, high pupil + locomotion


addpath(genpath('../pre_processing'));
addpath(genpath('../meta_data_processing'));
addpath(genpath('../../utils'));
datapath = 'X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData';
spike2path = 'X:\CardinLab\Antara\AnalyzedData\';
animals_db = get_animals_meta_data_by_csv;
procdatapath = 'X:\Hadas\Meso-imaging\crispr\meso_results\ProcessingDirectory';
%% Detects 3 arousal states per animal
for k=1:length(animals_db.folder_list)
    extract_sustained_state(datapath, procdatapath, spike2path, animals_db.folder_list{k});
end
end



function extract_sustained_state(datapath, procdatapath, spike2path, animalpath)
params.fspupilcam=30; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=5;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=5;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=3;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal)
params.minSitDuration=3;%minimum sit duration during quiescnece state

mkNewDir(fullfile(procdatapath, animalpath));
disp(animalpath);
spike2path = fullfile(spike2path, animalpath);
pupilfile = dir(fullfile(spike2path, '*_proc.mat'));
if  isempty(pupilfile)
    disp('no pupil file');
    return;
end
if  exist(fullfile(procdatapath, animalpath, 'arousal_state_ITI_segemts.mat'), 'file')
    disp('file exists');
    return;
end
if ~exist(fullfile(datapath, animalpath, 'smrx_signals_v2.mat'), 'file')
    warning('No smrx_signals_v2 file');
    return;
end

%% for each mouse, load spont and airpuff folders and perform correlations
load(fullfile(datapath, animalpath, 'smrx_signals_v2.mat'),'channels_data',...
    'timing');
if ~isfield(timing, 'airpuffstart')
    if isfield(timing, 'stimstart') % if this is a vis session, save vis as airpuff so we'll know to avoid these segmenents
        timing.airpuffstart = timing.stimstart;
        timing.airpuffend = timing.stimend;
    else
        timing.airpuffstart = [];
        timing.airpuffend = [];
    end
end
t_imaging = timing.bluestart;
dat=load(fullfile(pupilfile.folder, pupilfile.name));
if ~isfield(dat.proc, 'pupil')
    warning('No Pupil on face map');
    return;
end
pupil_Norm = dat.proc.pupil.area;
pupil_Norm = pupil_Norm - mean(pupil_Norm);
pupil_Norm = pupil_Norm/std(pupil_Norm);
pupil_time = timing.pupilcamstart;
if length(pupil_time) ~= length(pupil_Norm)
    warning('Pupil cam ticks are bad');
    return;
end
wheel_speed = channels_data.wheelspeed;
wheel_time = (1:length(wheel_speed))/params.fsspike2;

%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion s since locomotion
%onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)


[~, ~, ~, wheelOn, wheelOff] = get_wheel_on_off(raw_channels.data(:,5).', params.fsspike2, wheel_time, timing);

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
    return;
end
allEvtsStart=sort(timing.airpuffstart,'ascend')';
allEvtsEnd=sort(timing.airpuffend,'ascend')';

allEvtsPre=allEvtsStart-params.TimeSinceEvent;
allEvtsPost=allEvtsEnd+params.TimeSinceEvent;

newwheelOn=sort([wheelOn_t1,allEvtsPre,allEvtsPost],'ascend');
newwheelOff=sort([wheelOff_t1,allEvtsPre,allEvtsPost],'ascend');

Index=nan(1,length(newwheelOn));
for r=1:length(newwheelOn)
    tmp1=newwheelOn(r);
    tmp2=newwheelOff(r);
    
    tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
    if sum(tmp3)==0
        Index(r)=r;
    end
end
wheelOn_int=newwheelOn(Index(~isnan(Index)));
wheelOff_int=newwheelOff(Index(~isnan(Index)));

%makes sure the state is at least as long as the minimum run duration
idx1=find((wheelOff_int-wheelOn_int)>=(minRunDur));
wheelOn_int1=wheelOn_int(idx1);
wheelOff_int1=wheelOff_int(idx1);

%finalize the times to get sustained state only
wheelOn_final=wheelOn_int1+params.TimeSinceLocOn;
wheelOff_final=wheelOff_int1-params.TimeBeforeLocOff;

%% queiscence should be at least some criterion s long with some criterion s since locomotion offset
%and some criterion s before subsequent locomotion onset,excluding any events (airpuff/stim)
sitOn=[0;wheelOff(1:end-1)']; %use 0 as the first sit on time;
sitOff=wheelOn';%use wheelOn times as sit off times;

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

sitOn_t1=(tmpOn(:))'; sitOff_t1=(tmpOff(:))';

%find sit on off  times when airpuffs are not given

allEvtsStart=sort(timing.airpuffstart,'ascend')';
allEvtsEnd=sort(timing.airpuffend,'ascend')';

allEvtsPre=allEvtsStart-params.TimeSinceEvent;
allEvtsPost=allEvtsEnd+params.TimeSinceEvent;


newSitOn=sort([sitOn_t1,allEvtsPre,allEvtsPost],'ascend');
newSitOff=sort([sitOff_t1,allEvtsPre,allEvtsPost],'ascend');

Index=nan(1,length(newSitOn));
for r=1:length(newSitOn)
    tmp1=newSitOn(r);
    tmp2=newSitOff(r);
    
    tmp3= sum(tmp1>=allEvtsPre & tmp2<=allEvtsPost);
    if sum(tmp3)==0
        Index(r)=r;
    end
end
sitOn_int=newSitOn(Index(~isnan(Index)));
sitOff_int=newSitOff(Index(~isnan(Index)));


%makes sure state is at least as long as the minimum sit duration
idx1=find((sitOff_int-sitOn_int)>=(minSitDur));
sitOn_int1=sitOn_int(idx1);
sitOff_int1=sitOff_int(idx1);

%finalize the times to get sustained state only
sitOn_final=sitOn_int1+params.TimeSinceSitOn;
sitOff_final=sitOff_int1-params.TimeBeforeSitOff;
if 0
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
if ~isempty(allEvtsStart)
    for tt=1:length(allEvtsStart)
        plot([allEvtsStart(tt),allEvtsStart(tt)], ylimits, 'k');
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
    for tt=1:length(allEvtsStart)
        plot([allEvtsStart(tt),allEvtsStart(tt)], ylimits, 'k');
    end
end
plot([t_imaging(1),t_imaging(1)], ylimits, 'm');
plot([t_imaging(end),t_imaging(end)], ylimits, 'm');
if 0
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
    linkaxes([ax1, ax2, ax3,ax4,ax5,ax6],'x');
end

close all;
% Hadas/lav, change this once we get pupil data

segments_arousals.low_pup_q = [sitOn_final; sitOff_final]';
segments_arousals.high_pup_q = [];
segments_arousals.low_pup_q = [wheelOn_final; wheelOff_final]';

%     segments_arousals.low_pup_q = [Pupil_LowArousal_On_final_qu Pupil_LowArousal_Off_final_qu];
%     segments_arousals.high_pup_q = [Pupil_HighArousal_On_final_qu Pupil_HighArousal_Off_final_qu];
%     segments_arousals.high_pup_l = [Pupil_HighArousal_On_final_loc Pupil_HighArousal_Off_final_loc];
save(fullfile(procdatapath, animalpath, 'arousal_state_ITI_segemts.mat'),...
    'segments_arousals');



end