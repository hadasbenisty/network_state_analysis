function [timing, channels_data] = process_spike2_two_cams(cedpath, outputpath,data_smr_time_stamp_filename, fsspike2,channels,minRunDuration,minSitDuration,ITITime,pupilSR,tiffsPath)
%outputs timing in second and continuous data 
%load SON library
addpath(genpath('../pupil_proc'));
CEDS64LoadLib(cedpath)
channels_num=length(fieldnames(channels)); 
% parameters for extracting wheel motion
sCFG.sPARAM.dbWheelDiameterM=0.1524;%in m 
sCFG.sPARAM.db1WheelVRange=[0 5];%cut off range for signal
sCFG.sPARAM.dbWindowLenSec=2; % set the temoporal resolution of the analysis, this seems optimal so far
sCFG.sPARAM.blDoPlot=1; % indicate whether you want to plot or not
sCFG.sPARAM.blDoPlotDuration=1:500000; %plot duration in sample units
sCFG.sPARAM.winlengthsmoothing = 1;%smoothing wheel with 1 second window
%% load event channels 
data = process_spike2_smr2mat('', outputpath, data_smr_time_stamp_filename, channels_num);
[channels_data,wheelOn,wheelOff,h1] = get_channels_data_from_samples_with_wheelonoff(data, channels, sCFG,fsspike2);
mysave(h1, char(strcat(outputpath,'\wheel_on_off_times')), 'all');
timing.allwheelon=wheelOn(:); timing.allwheeloff=wheelOff(:);
threshold1=0.6;threshold2=0.05;
[timing.mesostart,timing.mesoend]=squaredetect(channels_data.mesoframe,.5);%camera start, end
timing.mesostart=timing.mesostart/fsspike2; timing.mesoend=timing.mesoend/fsspike2; 
[timing.bluestart,timing.blueend]=squaredetect(channels_data.blue,.5);%blue light start, end
timing.bluestart=timing.bluestart/fsspike2;timing.blueend=timing.blueend/fsspike2; 
[timing.uvstart,timing.uvend]=squaredetect(channels_data.uv,.5);%uv light start, end 
timing.uvstart=timing.uvstart/fsspike2; timing.uvend=timing.uvend/fsspike2; 
[airpuffstart,airpuffend]=squaredetect(channels_data.air_puff,.5);%airpuf/electrical stim start,end
airpuffstart=airpuffstart/fsspike2;airpuffend=airpuffend/fsspike2; 
if ~isempty(airpuffstart)
tmp=find(diff(airpuffstart)>5)+1; %make sure airpuffs/estim are spaced by at least 5s 
tmp1=airpuffstart(tmp);tmp2=airpuffend(tmp); 
airpuffstart=[airpuffstart(1); tmp1];
airpuffend=[airpuffend(1); tmp2];%this airpuff end time is not accurate however but we don't need this time 
end 
timing.allairpuffstart=airpuffstart(:); timing.allairpuffend=airpuffend(:); 
threshold1=0.6;%lav changed threshold 2 from 0.05 to 1 on dec 22 2020
threshold2=0.05;

[timing.pupilcamstart,timing.pupilcamend]= PupilDetectOnOff(channels_data.pupil,threshold1,threshold2,1/pupilSR,fsspike2);
if isempty(timing.mesostart)
timing.mesostart = timing.bluestart;
timing.mesoend = timing.blueend;
end
%% finetune locomotion and airpuff/electrical stimulation timing 
%only use wheel on with quiescence period of 10s and at least 5s of
%running and not occuring during other events (visual stim and airpuff)