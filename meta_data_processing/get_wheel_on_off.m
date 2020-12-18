function [wheelOn_int, wheelOff_int, wheelspeed, wheelOn, wheelOff] = get_wheel_on_off(subsetwheeldata, fsspike2, dataWheeltime, timing)

dataWheel.time{1} = dataWheeltime;
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
%             mysave(h1, fullfile(datapath, [animalNames{animal_i} num2str(days(day_i)) 'wheel_on_off']), 'fig');
wheelspeed=sCFG.sL0PPWR.db1SpeedMpS;
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

wheelOn_int=wheelOn_t1(idx3);
wheelOff_int=wheelOff_t1(idx3);