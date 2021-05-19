function pre_process_spont_animals_crispr
%%Code to process the ITI period for all animals with extraction
%of 3 arousal states: low pupil + quiescence, high pupil + quiescence, high pupil + locomotion
% or 2 states: quiescence, locomotion

addpath(genpath('..\network_state_analysis\pre_processing'));
addpath(genpath('..\network_state_analysis\meta_data_processing'));
addpath(genpath('..\..\network_state_analysis\utils'));
spike2path = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
animals_db = get_animals_meta_data_by_csv;
procdatapath = 'X:\Hadas\Meso-imaging\CRISPR\analysis_results';
mkNewDir(procdatapath);
doover=true;



for k=1:length(animals_db.folder_list)
    disp(k)
    isimaging_good =animals_db.toinclude_list(k) == find(strcmp(animals_db.toinclude_lut,'Good'));
    if isimaging_good
        extract_sustained_states_wheel_pupil_face(spike2path, procdatapath, animals_db.folder_list{k});
    end
    close all;
end

%
%  for k=1:length(animals_db.folder_list)
%     disp(k)
%     isimaging_good =animals_db.toinclude_list(k) == find(strcmp(animals_db.toinclude_lut,'Good'));
%     if isimaging_good
%         %% 2 arousal states per animal
%         resfile = fullfile(procdatapath, animals_db.folder_list{k}, 'arousal_2state_ITI_segemts.mat');
%         if ~isfile(resfile)|| doover
%             extract_sustained_2states(spike2path, procdatapath, animals_db.folder_list{k});
%         end
%         %% 3 arousal states per animal
%         resfile = fullfile(procdatapath, animals_db.folder_list{k}, 'arousal_3state_ITI_segemts.mat');
%         if ~isfile(resfile)|| doover
%             try
%             extract_sustained_3states(spike2path, procdatapath, animals_db.folder_list{k});
%             catch
%                 disp(['Failed to extract pupil for ' animals_db.folder_list{k}]);
%             end
%         end
%     end
%     close all;
% end
% 2 arousal states per animal
dffpath = 'X:\Hadas\Meso-imaging\CRISPR\traces_data';
% concatenateSpontPeriodsBy2States(doover, animals_db, dffpath, spike2path, procdatapath)
% 3 arousal states per animal
%concatenateSpontPeriodsByState_all(doover, animals_db, dffpath, spike2path, procdatapath)
concatenateSpontPeriodsByState(doover, animals_db, dffpath, spike2path, procdatapath)



disp('Pre prossing finished');
end
function extract_sustained_states_wheel_pupil_face(tracespath, procdatapath, animalpath)
params.fspupilcam=10; %pupil sampling rate
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=1;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=1;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=5;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=5;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=3;%for any state, minimum time since any event onset/offset
params.minRunDuration=5;% minimum run duration during locomotion state
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal)
params.minSitDuration=3;%minimum sit duration during quiescnece state

mkNewDir(fullfile(procdatapath, animalpath));
disp(animalpath);

if isfile(fullfile(procdatapath, animalpath, 'arousal_traces_states.mat'))
    disp('exist');
    return;
end
%% load spike2
if ~exist(fullfile(tracespath, animalpath, 'smrx_signals_v4.mat'), 'file')
    warning('No smrx_signals_v4 file');
    return;
end
load(fullfile(tracespath, animalpath, 'smrx_signals_v4.mat'),'channels_data',...
    'timing');
timing.events2ignorestart = [];
timing.events2ignoreend = [];
if isfield(timing, 'stimstart') % if this is a vis session
    timing.events2ignorestart = reshape(timing.stimstart, 1, []);
    timing.events2ignoreend = reshape(timing.stimend, 1, []);
end

if isfield(timing, 'airpuffstart')
    timing.events2ignorestart = cat(2, timing.events2ignorestart,reshape(timing.airpuffstart, 1, []));
    timing.events2ignoreend = cat(2, timing.events2ignoreend,reshape(timing.airpuffend, 1, []));
end

wheel_speed = channels_data.wheelspeed;
wheel_time = (1:length(wheel_speed))/params.fsspike2;

t_imaging = timing.bluestart;
%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion s since locomotion
%onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)


wheelOn=timing.wheelOn;
wheelOff = timing.wheelOff;

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
allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent,1,[]);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);

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
sitOn=([0;wheelOff(1:end-1)])'; %use 0 as the first sit on time;
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
tokeep = [];
for rr=1:length(sitOn_t1)
    st = findClosestDouble(wheel_time, sitOn_t1(rr)+0.1);
    en = findClosestDouble(wheel_time, sitOff_t1(rr)-0.1);
    if any(abs(wheel_speed(st:en))>speedThres)
        tokeep(rr) = false;
    else
        tokeep(rr) = true;
    end
end
sitOn_t1=sitOn_t1(tokeep==true);
sitOff_t1=sitOff_t1(tokeep==true);
%% this is Sweyta's code, could not understand it

% tmpOn=sitOn_t1; tmpOff=sitOff_t1;
% for rr=1:length(sitOn_t1)
%     currspeed = wheel_speed(wheel_time>sitOn_t1(rr)& wheel_time<sitOff_t1(rr));
%     highSpeedIdx=abs(currspeed)>speedThres;
%     if ~isempty(highSpeedIdx)
%         firstIdx=wheel_time(find(wheel_time>sitOn_t1(rr),1)+highSpeedIdx(1))-0.1;
%         lastIdx=wheel_time(find(wheel_time>sitOn_t1(rr),1)+highSpeedIdx(end))+0.1;
%         tmpOff(rr)=firstIdx;
%         tmpOff(end+1)=sitOff_t1(rr);
%         tmpOn(end+1)=lastIdx;
%     end
% end
% tmpOn=sort(tmpOn,'ascend'); tmpOff=sort(tmpOff,'ascend');
%
% sitOn_t1=(tmpOn(:))'; sitOff_t1=(tmpOff(:))';

%find sit on off  times when airpuffs are not given

allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent, 1, []);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);


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
segments_arousals.sit =  [sitOn_final' sitOff_final'];
segments_arousals.loc =  [wheelOn_final' wheelOff_final'];




%% pupil
ispupil = true;
pupilfile = dir(fullfile(tracespath,animalpath, 'pupil_clean.mat'));
if  isempty(pupilfile) || length(timing.pupilcamstart) < 1000
    disp('no pupil file');
    ispupil = false;
    pupil_Norm=nan;
else
    dat=load(fullfile(pupilfile.folder, pupilfile.name));
    if ~isfield(dat, 'areaii')
        warning('No Pupil on face map');
        return;
    end
    
    pupil = dat.areaii;
    pupil_Norm = pupil - mean(pupil);
    pupil_Norm = pupil_Norm/std(pupil_Norm);
    
    
    if size(pupil_Norm,1)==length(timing.pupilcamstart(1:end))
        timing.pupilcamstart=timing.pupilcamstart(1:end);
    elseif size(pupil_Norm,1)<length(timing.pupilcamstart(1:end))
        timing.pupilcamstart=timing.pupilcamstart(1:size(pupil_Norm,1));
        timing.pupilcamend=timing.pupilcamend(1:size(pupil_Norm,1));
    elseif size(pupil_Norm,1)>length(timing.pupilcamstart(1:end))
        pupil_Norm=pupil_Norm(1:length(timing.pupilcamstart),:);
        pupil=pupil(1:length(timing.pupilcamstart),:);
    end
    pupil_time=timing.pupilcamstart;
    %% do change point detection on pupil and face to get pupil high/low arousal or face high/low movement times during sustained quiescence state
    %get Z-thresholds based on pupil data during quiescence, when mouse
    %isn't moving and when aripuffs are not given
    b1DoPlot=1; blDoPlotDuration=1:800; smoothWin=1;
    pupilTime_Idx=cell(1,length(sitOn_int));
    for st=1:length(sitOn_int)
        pupilTime_Idx{st}=find(pupil_time>sitOn_int(st) & pupil_time <sitOff_int(st));
    end
    pupilTime_quiescence=cell2mat(pupilTime_Idx');
    pupil_qui_times=pupil_time(pupilTime_quiescence);
    pupil_quiescence=pupil_Norm(pupilTime_quiescence);
    zthres_High=quantile(pupil_quiescence,0.60);
    zthres_Low=quantile(pupil_quiescence,0.40);
    if all(isnan(zthres_High)) || all(isnan(zthres_Low))
        ispupil = false;
    else
        %get on and off timestamps for high and low arousal based on pupil
        [h1,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp ] =changepoints1(pupil_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
        [h2,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints1(-pupil_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
        
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
        segments_arousals.low_pupil = [Pupil_LowArousal_On_final_qu Pupil_LowArousal_Off_final_qu];
        segments_arousals.high_pupil = [Pupil_HighArousal_On_final_qu Pupil_HighArousal_Off_final_qu];
        
        
    end
end
%% face
isfacemap = true;
facemapfile = dir(fullfile(tracespath, animalpath,  'facemap.mat'));
if isempty(facemapfile) || ~ispupil
    warning('no facemap file')
    isfacemap = false;
else
    %get z thresholds based on face data during quiescence only
    dat=load(fullfile(tracespath, animalpath,  'facemap.mat'));
    face = dat.proc.motSVD{1}(:,1);
    face_Norm = zscore(dat.proc.motSVD{1}(:,1));
    b1DoPlot=1; blDoPlotDuration=1:4000; smoothWin=1;
    face_quiescence=face_Norm(pupilTime_quiescence);
    zthres_High=quantile(face_quiescence,0.60);
    zthres_Low=quantile(face_quiescence,0.40);
    
    %get on and off timestamps for high and low face movment
    [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints1(face_Norm', zthres_High,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,0);
    [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints1(-face_Norm', -zthres_Low,pupil_time,params.fspupilcam,smoothWin, b1DoPlot,blDoPlotDuration,1);
    
    
    if size(Face_LowArousal_OnTStamp,1)>size(Face_LowArousal_OffTStamp,1)
        Face_LowArousal_OnTStamp=Face_LowArousal_OnTStamp(1:size(Face_LowArousal_OffTStamp,1),:)
    else
    end
    
    
    % determine that pupil and face high/low arousal/movement time are at least minimum criterion seconds long
    idx1=find((Face_HighArousal_OffTStamp-Face_HighArousal_OnTStamp)>=params.minArousalDuration);
    Face_HighArousal_OnT_int=Face_HighArousal_OnTStamp(idx1);
    Face_HighArousal_OffT_int=Face_HighArousal_OffTStamp(idx1);
    idx2=find((Face_LowArousal_OffTStamp-Face_LowArousal_OnTStamp)>=params.minArousalDuration);
    Face_LowArousal_OnT_int=Face_LowArousal_OnTStamp(idx2);
    Face_LowArousal_OffT_int=Face_LowArousal_OffTStamp(idx2);
    
    %% Face on/off for q
    % get Face on times if both on and off times occur entirely during sustained queiscence states identified in the previous step
    toDelete=ones(1,length(Face_HighArousal_OnT_int));
    for rj=1:length(Face_HighArousal_OnT_int)
        tmp = find (Face_HighArousal_OnT_int(rj)>=sitOn_final & Face_HighArousal_OffT_int(rj)<=sitOff_final);
        toDelete(rj)=isempty(tmp);
    end
    Face_HighArousal_On_final_qu=Face_HighArousal_OnT_int(~toDelete);
    Face_HighArousal_Off_final_qu=Face_HighArousal_OffT_int(~toDelete);
    
    toDelete=ones(1,length(Face_LowArousal_OnT_int));
    for rj=1:length(Face_LowArousal_OnT_int)
        tmp = find (Face_LowArousal_OnT_int(rj)>=sitOn_final & Face_LowArousal_OffT_int(rj)<=sitOff_final);
        toDelete(rj)=isempty(tmp);
    end
    Face_LowArousal_On_final_qu=Face_LowArousal_OnT_int(~toDelete);
    Face_LowArousal_Off_final_qu=Face_LowArousal_OffT_int(~toDelete);
    %% pupil on/off for locomotion
    % get pupil on and face on times if both on and off times occur entirely
    %during sustained locomotion states identified in the previous step
    [Face_HighArousal_On_final_loc, Face_HighArousal_Off_final_loc,...
        Face_LowArousal_On_final_loc, Face_LowArousal_Off_final_loc] = ...
        getHighLowonOf4Locomotion(wheelOn_final, wheelOff_final, t_imaging, ...
        Face_HighArousal_OnT_int, Face_HighArousal_OffT_int, ...
        Face_LowArousal_OnT_int, Face_LowArousal_OffT_int);
    segments_arousals.low_face = [Face_LowArousal_On_final_qu Face_LowArousal_Off_final_qu];
    segments_arousals.high_face = [Face_HighArousal_On_final_qu Face_HighArousal_Off_final_qu];
    
    
end




%
% %% plot on and off times for sanity check
% %plot locomotion state on off times
% h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
% set(0,'CurrentFigure',h);ax1=subplot(6,1,1);
% plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('LocomotionState');
% for tt=1:length(wheelOn_final)
%     plot([wheelOn_final(tt),wheelOn_final(tt)], ylimits, 'g');
%     plot([wheelOff_final(tt),wheelOff_final(tt)], ylimits, 'r');
% end
% %plot airpuff timea andn imaging onset/offset
% if ~isempty(allEvtsStart)
%     for tt=1:length(allEvtsStart)
%         plot([allEvtsStart(tt),allEvtsStart(tt)], ylimits, 'k');
%     end
% end
% plot([t_imaging(1),t_imaging(1)], ylimits, 'm');
% plot([t_imaging(end),t_imaging(end)], ylimits, 'm');
%
% %plot quiescence state on off times
% set(0,'CurrentFigure',h);ax2=subplot(6,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
% for tt=1:length(sitOn_final)
%     plot([sitOn_final(tt),sitOn_final(tt)], ylimits, 'g');
%     plot([sitOff_final(tt),sitOff_final(tt)], ylimits, 'r');
% end
% %plot airpuff timea and imaging onset/offset times
% if exist('allEvts','var')
%     for tt=1:length(allEvtsStart)
%         plot([allEvtsStart(tt),allEvtsStart(tt)], ylimits, 'k');
%     end
% end
% plot([t_imaging(1),t_imaging(1)], ylimits, 'm');
% plot([t_imaging(end),t_imaging(end)], ylimits, 'm');
%
% %plot pupil and face state on off times
% set(0,'CurrentFigure',h);ax3=subplot(6,1,3);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilHighState+q');
% for tt=1:length(Pupil_HighArousal_On_final_qu)
%     plot([Pupil_HighArousal_On_final_qu(tt),Pupil_HighArousal_On_final_qu(tt)], ylimits, 'g');
%     plot([Pupil_HighArousal_Off_final_qu(tt),Pupil_HighArousal_Off_final_qu(tt)], ylimits, 'r');
% end
%
% set(0,'CurrentFigure',h);ax4=subplot(6,1,4);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilLowState+q');
% for tt=1:length(Pupil_LowArousal_On_final_qu)
%     plot([Pupil_LowArousal_On_final_qu(tt),Pupil_LowArousal_On_final_qu(tt)], ylimits, 'g');
%     plot([Pupil_LowArousal_Off_final_qu(tt),Pupil_LowArousal_Off_final_qu(tt)], ylimits, 'r');
% end
% %plot pupil and face state on off times
% set(0,'CurrentFigure',h);ax5=subplot(6,1,5);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilHighState+loc');
% for tt=1:length(Pupil_HighArousal_On_final_loc)
%     plot([Pupil_HighArousal_On_final_loc(tt),Pupil_HighArousal_On_final_loc(tt)], ylimits, 'g');
%     plot([Pupil_HighArousal_Off_final_loc(tt),Pupil_HighArousal_Off_final_loc(tt)], ylimits, 'r');
% end
%
% set(0,'CurrentFigure',h);ax6=subplot(6,1,6);plot(pupil_time,pupil_Norm);ylimits=ylim; hold on; title('PupilLowState+loc');
% for tt=1:length(Pupil_LowArousal_On_final_loc)
%     plot([Pupil_LowArousal_On_final_loc(tt),Pupil_LowArousal_On_final_loc(tt)], ylimits, 'g');
%     plot([Pupil_LowArousal_Off_final_loc(tt),Pupil_LowArousal_Off_final_loc(tt)], ylimits, 'r');
% end
% linkaxes([ax1, ax2, ax3,ax4,ax5,ax6],'x');
%
% mysave(gcf, fullfile(procdatapath, animalpath, 'arousal_sustained_state_traces'))



wheel_speed = interp1(wheel_time, wheel_speed, t_imaging);
if ispupil
    pupil = interp1(pupil_time, pupil, t_imaging);
else
    pupil=nan;
end
if isfacemap
    face = interp1(pupil_time, face, t_imaging);
else
    face=nan;
end

save(fullfile(procdatapath, animalpath, 'arousal_traces_states.mat'),...
    't_imaging','pupil','face', 'wheel_speed', 'segments_arousals');


if ~isfield(segments_arousals, 'low_face')||~isfield(segments_arousals, 'high_face')
    disp('no facemap fields')
elseif isempty(segments_arousals.low_face)||isempty(segments_arousals.high_face)
    disp('empty facemap')
else
end
    
figure;
subplot(3,1,1);plot(t_imaging, wheel_speed);
for ti=1:size(segments_arousals.loc,1)
    line([1 1]*segments_arousals.loc(ti,1), get(gca,'YLim'),'Color','g');
    line([1 1]*segments_arousals.loc(ti,2), get(gca,'YLim'),'Color','g');
end
for ti=1:size(segments_arousals.sit,1)
    line([1 1]*segments_arousals.sit(ti,1), get(gca,'YLim'),'Color','r');
    line([1 1]*segments_arousals.sit(ti,2), get(gca,'YLim'),'Color','r');
end
if ispupil
    subplot(3,1,2);plot(t_imaging, pupil);
    for ti=1:size(segments_arousals.high_pupil,1)
        line([1 1]*segments_arousals.high_pupil(ti,1), get(gca,'YLim'),'Color','g');
        line([1 1]*segments_arousals.high_pupil(ti,2), get(gca,'YLim'),'Color','g');
    end
    for ti=1:size(segments_arousals.low_pupil,1)
        line([1 1]*segments_arousals.low_pupil(ti,1), get(gca,'YLim'),'Color','r');
        line([1 1]*segments_arousals.low_pupil(ti,2), get(gca,'YLim'),'Color','r');
    end
end
if isfacemap
    subplot(3,1,3);plot(t_imaging, face);
    for ti=1:size(segments_arousals.high_face,1)
        line([1 1]*segments_arousals.high_face(ti,1), get(gca,'YLim'),'Color','g');
        line([1 1]*segments_arousals.high_face(ti,2), get(gca,'YLim'),'Color','g');
    end
    for ti=1:size(segments_arousals.low_face,1)
        line([1 1]*segments_arousals.low_face(ti,1), get(gca,'YLim'),'Color','r');
        line([1 1]*segments_arousals.low_face(ti,2), get(gca,'YLim'),'Color','r');
    end
end
% mysave(gcf, fullfile(procdatapath, animalpath, 'arousal_traces_states'));
% 
end
function concatenateSpontPeriodsBy2States(dover, animals_db, dffpath, spike2path, procdatapath)
[~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
regionLabel.Allen = allen_parcels.regionNum;
regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
statesnames = {'qui','loc'};
for ir=1:length(animals_db.animal_list)
    animal=animals_db.folder_list{ir};
    disp(animal);
    resfile = fullfile(procdatapath,  animal, 'spont_data_2states_dfff.mat');
    isimaging_good = animals_db.toinclude_list(ir)==find(strcmp(animals_db.toinclude_lut, 'Good'));
    
    if ~isimaging_good
        continue;
    end
    if exist(resfile, 'file')&&~dover
        continue;
    end
    spike2pth = fullfile(spike2path, animal);
    
    
    for si = 1:length(statesnames)
        eval([statesnames{si} '.Allen = [];']);
        %         eval([statesnames{si} '.Grid4 = [];']);
        eval([statesnames{si} '.t = [];']);
    end
    
    %     fsimaing=10;
    %     delay_filt=150;
    datafile_allen = fullfile(dffpath, animal,  'Ca_traces_spt_patch11_Allen_dfff.mat');
    %     datafile_grid4 = fullfile(dffpath, animal,'Ca_traces_spt_patch11_Grid4_dfff.mat');
    if ~exist(fullfile(spike2pth, 'smrx_signals_v4.mat'), 'file')
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v4.mat'),'timing');
    t_imaging = timing.bluestart;
    
    
    segmentfile = fullfile(procdatapath,animal,  'arousal_2state_ITI_segemts.mat');
    if exist(datafile_allen,'file')  && exist(segmentfile,'file')%
        load(segmentfile, 'segments_arousals');
        a=load(datafile_allen);
        %         g4=load(datafile_grid4);
        %         if length(t_imaging) > size(g4.parcels_time_trace,2)
        %             t_imaging=t_imaging(1:size(g4.parcels_time_trace,2));
        %         end
        
        standardinds=load('X:\Hadas\Meso-imaging\lan\xspsych\spt\xs_31_grid4_dfff.mat','par_inds');
        [parcels_names1, regionLabel1, finalindex1, ~, maskByAllen1, roiLabelsbyAllen1] = getAllenClusteringLabelsGrid(standardinds.par_inds, 4);
        %         [parcels_names.grid4, regionLabel.grid4, finalindex.grid4, ~, maskByAllen.grid4, roiLabelsbyAllen.grid4] = getAllenClusteringLabelsGrid(g4.par_inds, 4);
        %         if length(finalindex.grid4) == length(finalindex1)
        %             Xg4  =g4.parcels_time_trace(finalindex.grid4,:);
        %         else
        
        %             Xg4 = nan(length(finalindex1), size(g4.parcels_time_trace,2));
        %
        %             for hi = 1:length(finalindex1)
        %                 indsel = find(finalindex.grid4 == finalindex1(hi));
        %                 if ~isempty(indsel)
        %                     Xg4( (hi), :) = g4.parcels_time_trace(indsel, :);
        %                 else
        %                     disp(finalindex1(hi))
        %                 end
        %             end
        %         end
        %         parcels_names.grid4=parcels_names1;
        %         regionLabel.grid4=regionLabel1;
        %         finalindex.grid4=finalindex1;
        %         maskByAllen.grid4=maskByAllen1;
        %         roiLabelsbyAllen.grid4=roiLabelsbyAllen1;
        if length(t_imaging) > size(a.parcels_time_trace,2)
            t_imaging=t_imaging(1:size(a.parcels_time_trace,2));
        end
        
        Xa  =a.parcels_time_trace(finalindex.Allen,:);
        Xa = bsxfun(@minus, Xa, nanmean(Xa,2));
        Xa = bsxfun(@rdivide, Xa, nanstd(Xa,[],2));
        for si = 1:length(statesnames)
            if isfield(segments_arousals, statesnames{si})
                data = eval(statesnames{si});
                
                for seg_i = 1:size(segments_arousals.(statesnames{si}),1)
                    
                    seg = segments_arousals.(statesnames{si})(seg_i,:);
                    ind1 = findClosestDouble(t_imaging, seg(1));
                    ind2 = findClosestDouble(t_imaging, seg(2));
                    tt=ind1:ind2;
                    finalinds = tt(~isnan(sum(Xa(:, tt))));
                    data.Allen = cat(2, data.Allen, Xa(:, finalinds));
                    %                     data.Grid4 = cat(2, data.Grid4, Xg4(:, finalinds));
                    data.t=cat(1,data.t,t_imaging(ind1:ind2));
                end
            end
            eval([statesnames{si} '= data;']);
        end
        
        
        
        
    end
    
    
    save(resfile,'qui','loc');
    
end


end
function concatenateSpontPeriodsByState_all(dover, animals_db, dffpath, spike2path, procdatapath)
% [~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
% regionLabel.Allen = allen_parcels.regionNum;
% regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
statesnames = {'sit','loc','low_pupil','high_pupil','high_face','low_face'};
for ir=1:length(animals_db.animal_list)
    animal=animals_db.folder_list{ir};
    disp(animal);
    resfile = fullfile(procdatapath,  animal, 'con_states.mat');
    
    %     isimaging_good =isimaginggood(animals_db, ir);
    %     ispup_good =ispupilgood(animals_db, ir);
    %     if ~isimaging_good || ~ispup_good
    %         continue;
    %     end
    if exist(resfile, 'file')&&~dover
        continue;
    end
    spike2pth = fullfile(spike2path, animal);
    
    
    
    
    %     fsimaing=10;
    %     delay_filt=150;
    datafile_allen = fullfile(spike2path, animal,  'Ca_traces_spt_patch11_Allen_dfff.mat');
    %     datafile_grid4 = fullfile(dffpath, animal,'Ca_traces_spt_patch11_Grid4_dfff.mat');
    if ~exist(fullfile(spike2pth, 'smrx_signals_v4.mat'), 'file')
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v4.mat'),'timing');
    t_imaging = timing.bluestart;
    
    
    segmentfile = fullfile(procdatapath,animal,  'arousal_traces_states.mat');
    if ~exist(datafile_allen,'file')  || ~exist(segmentfile,'file')
        continue;
    end%%&& exist(datafile_grid4,'file')
    load(segmentfile, 'segments_arousals');
    a=load(datafile_allen);
    %         g4=load(datafile_grid4);
    %         if length(t_imaging) > size(g4.parcels_time_trace,2)
    %             t_imaging=t_imaging(1:size(g4.parcels_time_trace,2));
    %         end
    
    
    %         [parcels_names.grid4, regionLabel.grid4, finalindex.grid4, ~, maskByAllen.grid4, roiLabelsbyAllen.grid4] = getAllenClusteringLabelsGrid(g4.par_inds, 4);
    
    
    Xa  =a.parcels_time_trace(finalindex.Allen,:);
    %         Xg4  =g4.parcels_time_trace(finalindex.grid4,:);
    if length(t_imaging) > size(Xa,2)
        t_imaging=t_imaging(1:size(Xa,2));
    end
    Y = extract_segment(t_imaging, Xa, segments_arousals.sit);
    
    Xa = bsxfun(@minus, Xa, quantile(Y',.0050)');
    tind = randperm(size(Xa,2));
    Xa_shuffled = Xa(:,tind);
%     Xa = bsxfun(@rdivide, Xa, nanstd(Xa,[],2));
    res = get_data_by_state(statesnames, Xa, t_imaging, segments_arousals);
    names = fieldnames(res);
    
    for k=1:length(names)
       eval([names{k} ' =  res.(names{k});']); 
    end
    save(resfile,'sit',    'loc','low_pupil','high_pupil','low_face','high_face'); 
end



end




function concatenateSpontPeriodsByState(dover, animals_db, dffpath, spike2path, procdatapath)
% [~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
% regionLabel.Allen = allen_parcels.regionNum;
% regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% cd('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory');
statesnames = {'low_face','high_face','loc'};
for ir=1:length(animals_db.animal_list)
    animal=animals_db.folder_list{ir};
    disp(animal);
    resfile = fullfile(procdatapath,  animal, 'arousal_traces_states.mat');
    
    isimaging_good =isimaginggood(animals_db, ir);
    ispup_good =ispupilgood(animals_db, ir);
    if ~isimaging_good || ~ispup_good
        continue;
    end
    if exist(resfile, 'file')&&~dover
        continue;
    end
    spike2pth = fullfile(spike2path, animal);
    
    
    for si = 1:length(statesnames)
        eval([statesnames{si} '.Allen = [];']);
        %         eval([statesnames{si} '.Grid4 = [];']);
        eval([statesnames{si} '.t = [];']);
    end
    
    %     fsimaing=10;
    %     delay_filt=150;
    datafile_allen = fullfile(dffpath, animal,  'Ca_traces_spt_patch11_Allen_dfff.mat');
    %     datafile_grid4 = fullfile(dffpath, animal,'Ca_traces_spt_patch11_Grid4_dfff.mat');
    if ~exist(fullfile(spike2pth, 'smrx_signals_v4.mat'), 'file')
        continue;
    end
    load(fullfile(spike2pth, 'smrx_signals_v4.mat'),'timing');
    t_imaging = timing.bluestart;
    
    
    segmentfile = fullfile(procdatapath,animal,  'arousal_3state_ITI_segemts.mat');
    if exist(datafile_allen,'file')  && exist(segmentfile,'file')%%&& exist(datafile_grid4,'file')
        load(segmentfile, 'segments_arousals');
        a=load(datafile_allen);
        %         g4=load(datafile_grid4);
        %         if length(t_imaging) > size(g4.parcels_time_trace,2)
        %             t_imaging=t_imaging(1:size(g4.parcels_time_trace,2));
        %         end
        
        
        %         [parcels_names.grid4, regionLabel.grid4, finalindex.grid4, ~, maskByAllen.grid4, roiLabelsbyAllen.grid4] = getAllenClusteringLabelsGrid(g4.par_inds, 4);
        
        
        Xa  =a.parcels_time_trace(finalindex.Allen,:);
        %         Xg4  =g4.parcels_time_trace(finalindex.grid4,:);
        if length(t_imaging) > size(Xa,2)
            t_imaging=t_imaging(1:size(Xa,2));
        end
        Xa = bsxfun(@minus, Xa, nanmean(Xa,2));
        Xa = bsxfun(@rdivide, Xa, nanstd(Xa,[],2));
        
        for si = 1:length(statesnames)
            if isfield(segments_arousals, statesnames{si})
                data = eval(statesnames{si});
                
                for seg_i = 1:size(segments_arousals.(statesnames{si}),1)
                    
                    seg = segments_arousals.(statesnames{si})(seg_i,:);
                    ind1 = findClosestDouble(t_imaging, seg(1));
                    ind2 = findClosestDouble(t_imaging, seg(2));
                    tt=ind1:ind2;
                    finalinds = tt(~isnan(sum(Xa(:, tt))));
                    data.Allen = cat(2, data.Allen, Xa(:, finalinds));
                    %                     data.Grid4 = cat(2, data.Grid4, Xg4(:, finalinds));
                    data.t=cat(1,data.t,t_imaging(ind1:ind2));
                end
            end
            eval([statesnames{si} '= data;']);
        end
        
        
        
    end
    
    
    save(resfile,'low_face', 'high_face','loc');
    
end


end
function extract_sustained_2states(spike2path, procdatapath, animalpath)
params.fsspike2=5000;% spike2 sampling rate
params.TimeSinceLocOn=3;%for locomotion state, minimum time since locomotion onset
params.TimeBeforeLocOff=3;%for locomotion state, minimum time before locomotion offset
params.TimeSinceSitOn=5;%for quiescence state, minimum time since quiescence onset
params.TimeBeforeSitOff=5;%for quiescence state, minimum time before quiescence offset
params.TimeSinceEvent=3;%for any state, minimum time since any event onset/offset
params.minRunDuration=3;% minimum run duration during locomotion state
params.minArousalDuration=2; %minimum face/pupil arousal state (high or low arousal)
params.minSitDuration=3;%minimum sit duration during quiescnece state

mkNewDir(fullfile(procdatapath, animalpath));
disp(animalpath);

% if exist(fullfile(procdatapath, animalpath, 'arousal_2state_ITI_segemts.mat'), 'file')
%     disp('file exists');
%     return;
% end
if ~exist(fullfile(spike2path, animalpath, 'smrx_signals_v4.mat'), 'file')
    warning('No smrx_signals_v4 file');
    return;
end

%% for each mouse, load spont and airpuff folders and perform correlations
load(fullfile(spike2path, animalpath, 'smrx_signals_v4.mat'),'channels_data',...
    'timing');
timing.events2ignorestart = [];
timing.events2ignoreend = [];
if isfield(timing, 'stimstart') % if this is a vis session
    timing.events2ignorestart = reshape(timing.stimstart, 1, []);
    timing.events2ignoreend = reshape(timing.stimend, 1, []);
end

if isfield(timing, 'airpuffstart')
    timing.events2ignorestart = cat(2, timing.events2ignorestart,reshape(timing.airpuffstart, 1, []));
    timing.events2ignoreend = cat(2, timing.events2ignoreend,reshape(timing.airpuffend, 1, []));
end

t_imaging = timing.bluestart;

wheel_speed = channels_data.wheelspeed;
wheel_time = (1:length(wheel_speed))/params.fsspike2;

%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion s since locomotion
%onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)


wheelOn=timing.wheelOn;
wheelOff = timing.wheelOff;

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
allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent,1,[]);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);

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
sitOn=([0;wheelOff(1:end-1)])'; %use 0 as the first sit on time;
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
    highSpeedIdx=find(abs(currspeed)>speedThres);
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

allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent, 1, []);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);


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


%% plot on and off times for sanity check
%plot locomotion state on off times
h=figure; %initialize a figure that will contain subplots of behavior with on/off times markers for identified states
set(0,'CurrentFigure',h);ax1=subplot(2,1,1);
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
set(0,'CurrentFigure',h);ax2=subplot(2,1,2);plot(wheel_time,wheel_speed);ylimits=ylim; hold on; title('QuiescenceState');
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

linkaxes([ax1, ax2],'x');
% mysave(gcf, fullfile(procdatapath, animalpath, 'arousal_2state_traces'))



segments_arousals.loc = [wheelOn_final(:) wheelOff_final(:)];
segments_arousals.qui = [sitOn_final(:) sitOff_final(:)];
wheel_speedN = zscore(wheel_speed);
% extract pupil traces by state and concatenate
statesnames = fieldnames(segments_arousals);
for statei=1:length(statesnames)
    wheelvals.(statesnames{statei})=[];
    for k=1:size(segments_arousals.(statesnames{statei}),1)
        
        t1=findClosestDouble(segments_arousals.(statesnames{statei})(k,1), wheel_time);
        t2=findClosestDouble(segments_arousals.(statesnames{statei})(k,2), wheel_time);
        wheelvals.(statesnames{statei}) = cat(1, wheelvals.(statesnames{statei}), ...
            wheel_speedN(t1:t2)');
        
    end
end
save(fullfile(procdatapath, animalpath, 'arousal_2state_ITI_segemts.mat'),...
    'segments_arousals','wheelvals');




end

function extract_sustained_3states(spike2path, procdatapath, animalpath)
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

mkNewDir(fullfile(procdatapath, animalpath));
disp(animalpath);

pupilfile = dir(fullfile(spike2path,animalpath, 'pupil_clean.mat'));
if  isempty(pupilfile)
    disp('no pupil file');
    return;
end
if  exist(fullfile(procdatapath, animalpath, 'arousal_state_ITI_segemts.mat'), 'file')
    disp('file exists');
    return;
end
if ~exist(fullfile(spike2path, animalpath, 'smrx_signals_v4.mat'), 'file')
    warning('No smrx_signals_v4 file');
    return;
end

%% for each mouse, load spont and airpuff folders and perform correlations
load(fullfile(spike2path, animalpath, 'smrx_signals_v4.mat'),'channels_data',...
    'timing');
timing.events2ignorestart = [];
timing.events2ignoreend = [];
if isfield(timing, 'stimstart') % if this is a vis session
    timing.events2ignorestart = reshape(timing.stimstart, 1, []);
    timing.events2ignoreend = reshape(timing.stimend, 1, []);
end

if isfield(timing, 'airpuffstart')
    timing.events2ignorestart = cat(2, timing.events2ignorestart,reshape(timing.airpuffstart, 1, []));
    timing.events2ignoreend = cat(2, timing.events2ignoreend,reshape(timing.airpuffend, 1, []));
end

t_imaging = timing.bluestart;
dat=load(fullfile(pupilfile.folder, pupilfile.name));


if ~isfield(dat, 'areaii')
    warning('No Pupil on face map');
    return;
end

pupil_Norm = dat.areaii;
pupil_Norm = pupil_Norm - mean(pupil_Norm);
pupil_Norm = pupil_Norm/std(pupil_Norm);

if size(pupil_Norm,1)==length(timing.pupilcamstart(1:end))
    timing.pupilcamstart=timing.pupilcamstart(1:end);
elseif size(pupil_Norm,1)<length(timing.pupilcamstart(1:end))
    timing.pupilcamstart=timing.pupilcamstart(1:size(pupil_Norm,1));
    timing.pupilcamend=timing.pupilcamend(1:size(pupil_Norm,1));
elseif size(pupil_Norm,1)>length(timing.pupilcamstart(1:end))
    pupil_Norm=pupil_Norm(1:length(timing.pupilcamstart),:);
end
pupil_time=timing.pupilcamstart;

% if length(pupil_time) ~= length(pupil_Norm)
%     warning('Pupil cam ticks are bad');
%     return;
% end
disp('pupil cam ticks are good')

wheel_speed = channels_data.wheelspeed;
wheel_time = (1:length(wheel_speed))/params.fsspike2;

%% get locomotion on/off and quiescence on off times
%locomotion periods should be at least some criterion s long with some criterion s since locomotion
%onset, some criterion s before locomotion offset, excluding any events (airpuff/stim)


wheelOn=timing.wheelOn;
wheelOff = timing.wheelOff;

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
allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent,1,[]);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);

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
sitOn=([0;wheelOff(1:end-1)])'; %use 0 as the first sit on time;
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

allEvtsStart=sort(timing.events2ignorestart,'ascend')';
allEvtsEnd=sort(timing.events2ignoreend,'ascend')';

allEvtsPre=reshape(allEvtsStart-params.TimeSinceEvent, 1, []);
allEvtsPost=reshape(allEvtsEnd+params.TimeSinceEvent, 1, []);


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

%% do change point detection on pupil and face to get pupil high/low arousal or face high/low movement times during sustained quiescence state
%get Z-thresholds based on pupil data during quiescence, when mouse
%isn't moving and when aripuffs are not given
b1DoPlot=1; blDoPlotDuration=1:800; smoothWin=1;
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
if 0%~isempty(face_Norm)&exist(fullfile(procdatapath, animalpath, '*proc.mat'), 'file')
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

mysave(gcf, fullfile(procdatapath, animalpath, 'arousal_3state_traces'))



segments_arousals.low_pup_q = [Pupil_LowArousal_On_final_qu Pupil_LowArousal_Off_final_qu];
segments_arousals.high_pup_q = [Pupil_HighArousal_On_final_qu Pupil_HighArousal_Off_final_qu];
segments_arousals.high_pup_l = [Pupil_HighArousal_On_final_loc Pupil_HighArousal_Off_final_loc];
wheel_speedNorm = zscore(wheel_speed);
% extract pupil traces by state and concatenate
statesnames = fieldnames(segments_arousals);
for statei=1:length(statesnames)
    pupvals.(statesnames{statei})=[];
    wheelvals.(statesnames{statei})=[];
    for k=1:size(segments_arousals.(statesnames{statei}),1)
        t1=findClosestDouble(segments_arousals.(statesnames{statei})(k,1), pupil_time);
        t2=findClosestDouble(segments_arousals.(statesnames{statei})(k,2), pupil_time);
        pupvals.(statesnames{statei}) = cat(1, pupvals.(statesnames{statei}), ...
            pupil_Norm(t1:t2));
        
        t1=findClosestDouble(segments_arousals.(statesnames{statei})(k,1), wheel_time);
        t2=findClosestDouble(segments_arousals.(statesnames{statei})(k,2), wheel_time);
        wheelvals.(statesnames{statei}) = cat(1, wheelvals.(statesnames{statei}), ...
            wheel_speedNorm(t1:t2)');
        
    end
end

save(fullfile(procdatapath, animalpath, 'arousal_3state_ITI_segemts.mat'),...
    'segments_arousals','pupvals','wheelvals');
wheelN = interp1(wheel_time, wheel_speedNorm, pupil_time);
save(fullfile(procdatapath, animalpath, 'arousal_3state_traces.mat'),...
    'wheelN','pupil_time','pupil_Norm');


end