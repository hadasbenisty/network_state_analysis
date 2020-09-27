function main_organize_by_ITI_global_parcellation_lan(animalNames, days)
addpath(genpath('../parcellation/'));
addpath(genpath('../correlation'));
addpath(genpath('../LSSC-higley-master\LSSC-higley-master'));
addpath(genpath('../../../utils/Questionnaire/'));
addpath(genpath('../pre_processing_scripts'));
% load
[~, allen_parcels, maskAllen, maskAllenFrontal, allen_ROI_list] = getParcellsByLansAllansAtlass;


fltstr = 'pixelwise';
fltstr = 'spt';
%
%  animalNames = { 'xs'     };%'xu'   'xs'    };
%  days=[15]
fsspike2=5e3;
%initialize output and input paths
outfigs = 'X:\Hadas\Meso-imaging\lan';
spike2pth0 = 'X:\Hadas\Meso-imaging\lan\spike2data';
addpath(genpath('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master'));
dbstop if error;


data_smr_path = 'X:\Lan\Meso-imaging\';
cedpath = 'X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\pre_processing_scripts\utils\CEDS64ML';

for animal_i =1:length(animalNames)
    switch animalNames{animal_i}
        case {'xu','xv','xt','xs'}
            fsimaing=33;
            delay_filt = 500;
        otherwise
            fsimaing=10;
            delay_filt=150;
    end
%     par_gal_mask = load(['X:\Hadas\Meso-imaging\lan\' animalNames{animal_i} 'psych\' fltstr '/gal/' animalNames{animal_i} '_finalparcellation_Gal']);
    for day_i = 1:length(days)
        disp([animalNames{animal_i}  ' '  num2str(days(day_i))]);
        %load data 
        datapath = ['X:\Hadas\Meso-imaging\lan\' animalNames{animal_i} 'psych\' fltstr '\'];
        spike2pth = fullfile(spike2pth0, animalNames{animal_i});
        
        if exist(fullfile(datapath, [animalNames{animal_i} num2str(day_i) 'imaging_time_traces_global_ITI.mat']), 'file')
            continue;
        end
        % spike2
        if ~exist(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']),'file')
            filename = fullfile(data_smr_path, animalNames{animal_i}, [animalNames{animal_i} '_D' num2str(days(day_i))]);
            if ~exist([filename '.smrx'],'file')
                disp('No smrx');
                continue;
            end
            [timing, channels_data] = process_spike2(cedpath, '',filename, fsspike2);
            t_spike2 = linspace(0, length(channels_data.startsig)-1,length(channels_data.startsig))/fsspike2;
            t_imaging = timing.mesostart/fsspike2;
            
            save(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']), 't_spike2', 'channels_data', 't_imaging',...
                'timing', 't_spike2');
        else
            load(fullfile(spike2pth, ['spike2data',animalNames{animal_i} num2str(days(day_i)) '.mat']),'channels_data',...
                'timing', 't_imaging');
        end
        before_win = 1;
        after_win = 3;
        
        parfile = fullfile(datapath, [animalNames{animal_i} '_' num2str(days(day_i)) '_allen.mat']);
        if ~exist(parfile, 'file')
            disp('no allen');
            continue;
        end
        pardataAllan = load(parfile);
        
        %         parfileGal = fullfile(datapath,'gal',[animalNames{animal_i} '_' num2str(days(day_i)) ...
        %             '_global.mat']);
        %
        %
        %         if ~exist(parfileGal, 'file')
        %
        %             disp('no gal');
        %             continue;
        %         end
        %
        %         pardata_gal = load(parfileGal);
        %         [roiLabelsbyAllen.Gal, maskByAllen.Gal, regionLabel.Gal, isLeftLabel.Gal] = getAllenClusteringLabelsGlobal('', parfileGal, allen_parcels, maskAllenFrontal, par_gal_mask.final_par_mask);
        %
        
        
        if fsimaing < 33
            t_imaging=t_imaging(1:2:end);
        end
        t_imaging = t_imaging(1:end-delay_filt);
        
        if length(t_imaging) > size(pardataAllan.parcels_time_trace,2)
            t_imaging=t_imaging(1:size(pardataAllan.parcels_time_trace,2));
        end
        X_Allen =pardataAllan.parcels_time_trace(:, 1:length(t_imaging));
        %         X_Gal =pardata_gal.parcels_time_trace(:, 1:length(t_imaging));
        
        stim_timestamps = timing.stimstart/fsspike2;
        stim_timestamps=stim_timestamps(1:75);
        
        %         [~, imaging_time_traces.Gal] = time_trace2ITI(X_Gal, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        
        %get ITI periods between trials
        [imaging_time_traces.t, imaging_time_traces.Allen] = time_trace2ITI(X_Allen, t_imaging, stim_timestamps, before_win, after_win, fsimaing);
        figure;
        %plot timestamps of the first few seconds-> should contain up to
        %1-3 jumps as a sanity check the trial periods are being cut out
        plot(imaging_time_traces.t(1:131*fsimaing).');
        mysave(gcf, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'timestamps_imaging']), 'fig');
        close all;
        %Load region labels
        parcelsLabels.Allen = allen_parcels.regionNum;
        roiLabelsbyAllen.Allen = 1:length(allen_parcels.names);
        maskByAllen.Allen = maskAllen;
        regionLabel.Allen = allen_parcels.regionNum;
        isLeftLabel.Allen = repmat(1:2, 1, 28);
        regionLabel.nameslegend = {'Rest','Visual','Parietal','Temp','Aud','R-S','S-S','Motor'};
        
        save(fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'imaging_time_traces_global_ITI.mat']), 'imaging_time_traces' ,...
            'roiLabelsbyAllen', 'maskByAllen', 'regionLabel', 'isLeftLabel');
        clear parcelsLabels;
        clear trialslabels;
        t_spike2 = [1:length(channels_data.wheel)]';
        [ond,ofd] = squaredetect(channels_data.diode,.5);
        %get wheelspeed in the ITI period. use ond because t_spike2 is in
        %indices, not time and so is ond from squaredetect
        [~, spike2time_traces.wheelspeed] = time_trace2ITI(channels_data.wheelspeed.', t_spike2, ond, before_win, after_win, fsspike2);
        [spike2time_traces.t, spike2time_traces.wheel] = time_trace2ITI(channels_data.wheel.', t_spike2, ond, before_win, after_win, fsspike2);
        spike2time_traces.t=spike2time_traces.t/fsspike2;
        %spike2time_traces.t contains jumps where there are trial periods
        %and is the time vector for the wheelspeed in the ITI periods
        %% added by lav for wheel on and off timings
        raw_spike2 = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '.mat']);
        if ~exist(raw_spike2, 'file')
            disp('no raw_spike2');
            continue;
        end
        raw_channels = load(raw_spike2,  'data','timestamp'); %load raw wheel data from lan because we need the rotations for 
        %wheelspeed changepoints using the function below
        if length(raw_channels.timestamp)==length(t_spike2)&&size(raw_channels.data,1)==size(channels_data.wheel.',2) %ensure that the size of the raw channels data is the same as spike2 processed from hadas (so we can align them)
            [dataWheel.time{1}, subsetwheeldata] = time_trace2ITI(raw_channels.data(:,5).', raw_channels.timestamp, stim_timestamps, before_win, after_win, fsspike2);
            %wheelspeed in ITI periods
            %8/08/2020 changed spike2time_traces.t to raw_channels.timestamp and ond to
            %stim_timestamps. raw_channels.timestamp is effectively almost
            %the same as t_spike2 because previously it was incorrect. This
            %change was implemented for the data shown on 08/14/2020
            dataWheel.fsample=fsspike2;
            sCFG.sPARAM.dbWheelDiameterM=0.1524;%in m
            sCFG.sPARAM.db1WheelVRange=[0 5];%cut off range for signal
            sCFG.sPARAM.dbWindowLenSec=2; % set the temoporal resolution of the analysis, this seems optimal so far
            sCFG.sPARAM.blDoPlot=1; % indicate whether you want to plot or not 
            sCFG.sPARAM.blDoPlotDuration=1:500000; %plot duration in sample units
            %wheelData=data(:,channels.WHEEL);
            wheelData=fillmissing(subsetwheeldata,'nearest');
            dataWheel.trial{1}=(wheelData(:))';
            %dataWheel.time{1}=1/dataWheel.fsample:1/dataWheel.fsample:((size(dataWheel.trial{1},2))/dataWheel.fsample);
            hold off;
            [h1,sCFG] = wheel_changepoints(sCFG,dataWheel);
            mysave(h1, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'wheel_on_off']), 'fig');
            channels_data.wheelspeed=sCFG.sL0PPWR.db1SpeedMpS;
            wheelOn=sCFG.sL1DWCP.db1WheelOnTStamp;
            wheelOff=sCFG.sL1DWCP.db1WheelOffTStamp;
            
            minRunDuration=11;
            minSitDuration=15;
            %finetune locomotion parameters
            %only use wheel on with quiescence period of 15s and running
            %periods of 11 seconds because for running we cut out the first
            %3 and last 3 seconds and for sitting we cut out the first 5
            %and last 5 seconds
            idx=(wheelOn>(timing.mesostart(1)./fsspike2)+minSitDuration);
            wheelOn_t1=wheelOn(idx);
            wheelOff_t1=wheelOff(idx);
            %some indexing on the wheel on time stamps
            idx1=find((wheelOff_t1-wheelOn_t1)>=minRunDuration);
            newWheelOnTimes=wheelOn_t1(2:end);
            newWheelOffTimes=wheelOff_t1(1:end-1);
            idx2=(find((newWheelOnTimes-newWheelOffTimes)>=minSitDuration))+1;
            idx3=intersect(idx1,idx2);
            
            wheelOn_int=wheelOn_t1(idx3,1);
            wheelOff_int=wheelOff_t1(idx3,1);
            if length(wheelOn_int(:))>0
                tonotDelete_wheelon=ones(length(wheelOn_int(:)),1);tonotDelete_wheeloff=ones(length(wheelOff_int(:)),1);
                for rj=1:length(wheelOn_int)
                    %delete periods of onset and offset during stimulus time windows (sanity check)
                    tmp2=find(wheelOn_int(rj)>=stim_timestamps & wheelOff_int(rj)<=stim_timestamps+2*fsimaing);
                    tonotDelete_wheelon(rj)=isempty(tmp2); 
                end
                for rj=1:length(wheelOff_int)-1
                    %delete periods of onset and offset during stimulus time windows (sanity check)
                    tmp=find(wheelOff_int((rj))>=stim_timestamps & wheelOn_int((rj+1))<=stim_timestamps+2*fsimaing);
                    tonotDelete_wheeloff(rj)=isempty(tmp); 
                end
                finalwheelon=wheelOn_int(find(tonotDelete_wheelon==1));
                finalwheeloff=wheelOff_int(find(tonotDelete_wheeloff==1));
                %finalwheelOn=wheelOn_int(wheelOn_int>(dataWheel.time{1}(1)+ minSitDuration) & wheelOn_int<(dataWheel.time{1}(end)-minRunDuration));
                spike2time_traces.allwheelon=finalwheelon(:); spike2time_traces.allwheeloff=finalwheeloff(:);
                [running_time_traces.t_wheelon, running_time_traces.locomotionperiods] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allwheelon, fsimaing,8);
                %remove first 3 seconds and first 5 seconds (the inputs to the function
                %time_trace2ITI_running) ensures there are 5 seconds after
                %this deletion and only 5 seconds
                running_time_traces.t_wheelon(1:3*fsimaing,:)=[];
                running_time_traces.locomotionperiods(:,1:3*fsimaing,:)=[];
                [running_time_traces.t_wheeloff, running_time_traces.quiescenceperiods] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allwheeloff, fsimaing,10);
                running_time_traces.t_wheeloff(1:5*fsimaing,:)=[];
                running_time_traces.quiescenceperiods(:,1:5*fsimaing,:)=[];
                %remove any running periods that contain a temporal jump
                %due to cutting out trial periods
                toremove_locomotionperiods=[];
                for num=1:size(running_time_traces.t_wheelon,2)
                    tempvec=running_time_traces.t_wheelon(:,num);
                    if max(diff(tempvec)) >0.9
                        toremove_locomotionperiods(num)=1;
                    else
                        toremove_locomotionperiods(num)=0;
                    end
                end
                running_time_traces.t_wheelon(:,find(toremove_locomotionperiods==1))=[];
                running_time_traces.locomotionperiods(:,:,find(toremove_locomotionperiods==1))=[];
                
                toremove_quieperiods=[];
                for num=1:size(running_time_traces.t_wheeloff,2)
                    tempvec=running_time_traces.t_wheeloff(:,num);
                    if max(diff(tempvec)) >0.9
                        toremove_quieperiods(num)=1;
                    else
                        toremove_quieperiods(num)=0;
                    end
                end
                running_time_traces.t_wheeloff(:,find(toremove_quieperiods==1))=[];
                running_time_traces.quiescenceperiods(:,:,find(toremove_quieperiods==1))=[];
                
            else
            end
        else
            disp('rotation channel data not the same size as other channeldata');
        end
        %load in pupil data
        pupil_data = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_pupil_clean.mat']);
        if ~exist(pupil_data, 'file')
            disp('no pupil_data');
            continue;
        else
            pupilarea = load(pupil_data,  'areaii');
            pupil_Norm=pupilarea.areaii; %get pupil camera timepoints
            [timing.pupilstart,timing.pupilend]=squaredetect(channels_data.pupil_frame,.5);
            if size(pupil_Norm,1)==length(timing.pupilstart(1:end))
                timing.pupilstart=timing.pupilstart(1:end);
            elseif size(pupil_Norm,1)<length(timing.pupilstart(1:end))
                timing.pupilstart=timing.pupilstart(1:size(pupil_Norm,1));
                timing.pupilend=timing.pupilend(1:size(pupil_Norm,1));
            elseif size(pupil_Norm,1)>length(timing.pupilstart(1:end))
                pupil_Norm=pupil_Norm(1:length(timing.pupilstart),:);
            end
            pupil_time=timing.pupilstart/fsspike2;
            %%added by lav 08/06/2020 for extracting pupil and facemap
            %% do change point detection on pupil and face to get pupil on/off or whisker on/off times
            %from sweytas code
            b1DoPlot=1; blDoPlotDuration=100:4000; smoothWin=1;
            %get z thresholds based on pupil data during quiescence only
            pupilTime_qui=cell(1,length(wheelOff)-1);
            for st=1:length(wheelOff)-1
                pupilTime_qui{st}=find(pupil_time>wheelOff(st) & pupil_time <wheelOn(st+1));
            end
            pupilTime_quiescence=cell2mat(pupilTime_qui');
            running_time_traces.pupilTime_quiescence=pupilTime_quiescence;
            %running_time_traces.pupil_qui_times=pupil_time(pupilTime_quiescence);
            pupil_quiescence=pupil_Norm(pupilTime_quiescence); %subset pupil by the quiescence timestamps.
            
            running_time_traces.zthres_High=quantile(pupil_quiescence,0.60);
            running_time_traces.zthres_Low=quantile(pupil_quiescence,0.40); %get threshold
            
            %get on and off timestamps for high and low arousal based on pupil
            %pup_framerate=findClosestDouble(pupil_time,pupil_time(1)+1);
            pup_framerate=30;
            [h1,Pupil_HighArousal_OnTStamp,Pupil_HighArousal_OffTStamp] =changepoints(pupil_Norm', running_time_traces.zthres_High,pupil_time,pup_framerate,smoothWin, b1DoPlot,blDoPlotDuration);
            [h2,Pupil_LowArousal_OnTStamp,Pupil_LowArousal_OffTStamp ] =changepoints(-pupil_Norm', -running_time_traces.zthres_Low,pupil_time,pup_framerate,smoothWin, b1DoPlot,blDoPlotDuration);
            mysave(h1, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'pup_high']), 'fig');
            mysave(h2, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'pup_low']), 'fig');
            
            
            if length(Pupil_HighArousal_OnTStamp(:))>0&length(Pupil_LowArousal_OnTStamp(:))>0
                %remove locomotion and trial contaminated periods
                tonotDelete=ones(length(Pupil_HighArousal_OnTStamp(:)),1);
                for rj=1:min(length(Pupil_HighArousal_OnTStamp),length(Pupil_HighArousal_OffTStamp))
                    tmp = find (Pupil_HighArousal_OnTStamp(rj)>=wheelOn & Pupil_HighArousal_OffTStamp(rj)<=wheelOff); %find pupil onset and offset within a running periods 
                    tmp2=find (Pupil_HighArousal_OnTStamp(rj)>=stim_timestamps & Pupil_HighArousal_OffTStamp(rj)<=stim_timestamps+2*fsimaing);
                    tonotDelete(rj)=isempty(tmp)&isempty(tmp2); %if its empty (no running period pupil is found in), mark as 1 so its not deletted
                end
                Pupil_HighArousal_On_final=Pupil_HighArousal_OnTStamp(find(tonotDelete==1));
                Pupil_HighArousal_Off_final=Pupil_HighArousal_OffTStamp(find(tonotDelete==1));
                
                tonotDelete=ones(length(Pupil_LowArousal_OnTStamp),1);
                for rj=1:min(length(Pupil_LowArousal_OnTStamp),length(Pupil_LowArousal_OffTStamp ))
                    tmp = find (Pupil_LowArousal_OnTStamp(rj)>=wheelOn & Pupil_LowArousal_OffTStamp (rj)<=wheelOff);
                    tmp2=find (Pupil_LowArousal_OnTStamp(rj)>=stim_timestamps & Pupil_LowArousal_OffTStamp (rj)<=stim_timestamps+2*fsimaing);
                    tonotDelete(rj)=isempty(tmp)&isempty(tmp2);
                end
                Pupil_LowArousal_On_final=Pupil_LowArousal_OnTStamp(find(tonotDelete==1));
                Pupil_LowArousal_Off_final=Pupil_LowArousal_OffTStamp(find(tonotDelete==1));
               
                spike2time_traces.allpuphighon=Pupil_HighArousal_On_final(:)/30; spike2time_traces.allpuplowon=Pupil_LowArousal_On_final(:)/30;
                
                [running_time_traces.t_puphigh_on, running_time_traces.puphigh_on] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allpuphighon, fsimaing,5);
                [running_time_traces.t_puplow_on, running_time_traces.puplow_on] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allpuplowon, fsimaing,5);
                
                %%remove state transition if it contains a jump point (ie where signal
                %%jumps in time due to removal of trial or locomotion period
                toremove_puplow=[];
                for num=1:size(running_time_traces.t_puplow_on,2)
                    tempvec=running_time_traces.t_puplow_on(:,num);
                    if max(diff(tempvec)) >0.9
                        toremove_puplow(num)=1;
                    else
                        toremove_puplow(num)=0;
                    end
                end
                running_time_traces.t_puplow_on(:,find(toremove_puplow==1))=[];
                running_time_traces.puplow_on(:,:,find(toremove_puplow==1))=[];
                toremove_puphigh=[];
                for num=1:size(running_time_traces.t_puphigh_on,2)
                    tempvec=running_time_traces.t_puphigh_on(:,num);
                    if max(diff(tempvec)) >0.75
                        toremove_puphigh(num)=1;
                    else
                        toremove_puphigh(num)=0;
                    end
                end
                running_time_traces.t_puphigh_on(:,find(toremove_puphigh==1))=[];
                running_time_traces.puphigh_on(:,:,find(toremove_puphigh==1))=[];
                %ensure there is no overlap with (1) trial time and (2) locomotion times
                save(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features', [animalNames{animal_i} num2str(days(day_i)) 'running_ITI.mat']),  'running_time_traces');
            else
            end
        end
        %load face data
        facemap_data = fullfile('X:\Hadas\Meso-imaging\lan\facemap\Current_Processing\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_proc.mat']);
        if ~exist(facemap_data, 'file')
            disp('no facemap_data');
            continue;
        else
            facemapinfo = load(facemap_data);
            face_Norm=facemapinfo.proc.motSVD{1,1}(:,1);
            %face_Norm=medfilt1(facemapinfo.proc.motSVD{1,1}(:,1),3);
            face_Norm=fillmissing(face_Norm,'nearest');
            
            %facemap time
            [timing.pupilstart,timing.pupilend]=squaredetect(channels_data.pupil_frame,.5);
            
            if size(pupil_Norm,1)==length(timing.pupilstart(2:end))
                timing.pupilstart=timing.pupilstart(2:end);
            elseif size(pupil_Norm,1)<length(timing.pupilstart(2:end))
                timing.pupilstart=timing.pupilstart(2:(size(pupil_Norm,1)+1));
                timing.pupilend=timing.pupilend(2:(size(pupil_Norm,1)+1));
            elseif size(pupil_Norm,1)>length(timing.pupilstart(1:end))
                pupil_Norm=pupil_Norm(2:(length(timing.pupilstart)+1),:);
            end
            pupil_time=timing.pupilstart/fsspike2;
            
            %get z thresholds based on face data during quiescence only
            b1DoPlot=1; blDoPlotDuration=100:4000; smoothWin=1;
            face_quiescence=face_Norm(pupilTime_quiescence);
            zthres_High=quantile(face_quiescence,0.60);
            zthres_Low=quantile(face_quiescence,0.40);
            
            %get on and off timestamps for high and low face movment
            [h3,Face_HighArousal_OnTStamp,Face_HighArousal_OffTStamp] =changepoints(face_Norm', zthres_High,pupil_time,pup_framerate,smoothWin, b1DoPlot,blDoPlotDuration);
            [h4,Face_LowArousal_OnTStamp,Face_LowArousal_OffTStamp] =changepoints(-face_Norm', -zthres_Low,pupil_time,pup_framerate,smoothWin, b1DoPlot,blDoPlotDuration);
            mysave(h3, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'face_high']), 'fig');
            mysave(h4, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'face_low']), 'fig');
            
            if length(Face_HighArousal_OnTStamp(:))>0&length(Face_LowArousal_OnTStamp(:))>0                
                %%clean out times contaminated by locomotion and within the trial
                 tonotDelete=ones(length(Face_HighArousal_OnTStamp(:)),1);
                for rj=1:length(Face_HighArousal_OnTStamp)
                    tmp = find(Face_HighArousal_OnTStamp(rj)>=wheelOn & Face_HighArousal_OffTStamp(rj)<=wheelOff); %find Face onset and offset within a running periods 
                    tmp2=find(Face_HighArousal_OnTStamp(rj)>=stim_timestamps & Face_HighArousal_OffTStampp(rj)<=stim_timestamps+2*fsimaing);
                    tonotDelete(rj)=isempty(tmp)&isempty(tmp2); %if its empty (no running period Face is found in), mark as 1 so its not deletted
                end
                Face_HighArousal_On_final=Face_HighArousal_OnTStamp(find(tonotDelete==1));
                Face_HighArousal_Off_final=Face_HighArousal_OffTStamp(find(tonotDelete==1));
                
                tonotDelete=ones(length(Face_LowArousal_OnTStamp),1);
                for rj=1:length(Face_LowArousal_OnTStamp)
                    tmp = find(Face_LowArousal_OnTStamp(rj)>=wheelOn & Face_LowArousal_OffTStamp(rj)<=wheelOff);
                    tmp2=find(Face_LowArousal_OnTStamp(rj)>=stim_timestamps &Face_LowArousal_OffTStamp(rj)<=stim_timestamps+2*fsimaing);
                    tonotDelete(rj)=isempty(tmp)&isempty(tmp2);
                end
                Face_LowArousal_On_final=Face_LowArousal_OnTStamp(find(tonotDelete==1));
                Face_LowArousal_Off_final=Face_LowArousal_OffTStamp(find(tonotDelete==1));
               
                spike2time_traces.allfacehighon=Face_HighArousal_On_final(:)/30; spike2time_traces.allfacelowon=Face_LowArousal_On_final(:)/30;
                
                [running_time_traces.t_facehigh_on, running_time_traces.facehigh_on] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allfacehighon, fsimaing,5);
                [running_time_traces.t_facelow_on, running_time_traces.facelow_on] = time_trace2ITI_running(imaging_time_traces.Allen, imaging_time_traces.t, spike2time_traces.allfacelowon, fsimaing,5);
                %remove facemap timestamps with temporal jumps
                toremove_facelow=[];
                for num=1:size(running_time_traces.t_facelow_on,2)
                    tempvec=running_time_traces.t_facelow_on(:,num);
                    if max(diff(tempvec)) >0.9
                        toremove_facelow(num)=1;
                    else
                        toremove_facelow(num)=0;
                    end
                end
                running_time_traces.t_facelow_on(:,find(toremove_facelow==1))=[];
                running_time_traces.facelow_on(:,:,find(toremove_facelow==1))=[];
                toremove_facehigh=[];
                for num=1:size(running_time_traces.t_facehigh_on,2)
                    tempvec=running_time_traces.t_facehigh_on(:,num);
                    if max(diff(tempvec)) >0.9
                        toremove_facehigh(num)=1;
                    else
                        toremove_facehigh(num)=0;
                    end
                end
                running_time_traces.t_facehigh_on(:,find(toremove_facehigh==1))=[];
                running_time_traces.facehigh_on(:,:,find(toremove_facehigh==1))=[];
                
                disp('processed facemap');
                save(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features', [animalNames{animal_i} num2str(days(day_i)) 'running_ITI.mat']),  'running_time_traces');
            else
            end
        end
        
        clear imaging_time_traces;
        clear running_time_traces;
        %%
        
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features'));
        save(fullfile('X:\Hadas\Meso-imaging\lan\', [animalNames{animal_i} 'psych'], 'spike2Features', [animalNames{animal_i} num2str(days(day_i)) 'spike2_ITI.mat']),  'spike2time_traces');
        clear spike2time_traces;
        
        
        % whisking
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_whisker_clean.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisk_binary');
        wsk.whisk_binary =wsk.whisk_binary(1:min(length(wsk.whisk_binary),length(t_imaging)));
        [~, wsk_time_trace.binary] = time_trace2ITI(wsk.whisk_binary', t_imaging, stim_timestamps, before_win, after_win, pup_framerate);
        
        wskile = fullfile('X:\Lan\Meso-imaging\', animalNames{animal_i}, [animalNames{animal_i} '_D' ...
            num2str(days(day_i))  '_whisker.mat']);
        if ~exist(wskile, 'file')
            disp('no wskile');
            continue;
        end
        wsk = load(wskile,  'whisker');
        wsk.whisker =wsk.whisker(1:min(length(wsk.whisker),length(t_imaging)));
        
        [wsk_time_trace.t, wsk_time_trace.whisker] = time_trace2ITI(wsk.whisker', t_imaging, stim_timestamps, before_win, after_win, pup_framerate);
        
        mkdir(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures'))
        save(fullfile('X:\Hadas\Meso-imaging\lan\',[animalNames{animal_i} 'psych'], 'videoFeatures', [animalNames{animal_i} num2str(days(day_i)) 'whisk_ITI.mat']), ...
            'wsk_time_trace');
        clear 'channels_data';
        clear 'timing';
        
        
    end
    
    
    
end
