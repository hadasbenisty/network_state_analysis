function hannahnew_mecp21p_demixalign(Session,location,color,sessid,animal_group,cohort)

%% main preprocessing script for 1P analysis 
% session inputs as folder names of final session, location inputs as 'B449' or 'F238', color inputs as 'blueuv' or 'RCaMP_AC'
%and sessid as 'vis' or 'spont' depending on whether it's a visual stim
%session or a spontaneous session with locomotion and airpuff

if strcmp(animal_group,'emx')
    animal_group=strrep(animal_group,'emx','Emx');
elseif strcmp(animal_group,'ctrl')
    animal_group=strrep(animal_group,'ctrl','Control');
elseif strcmp(animal_group,'lhx6')
    animal_group=strrep(animal_group,'lhx6','Lhx6');
end

Session_after=strcat(animal_group,'_',Session);
%% add path 
addpath(genpath('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master'));
%addpath(genpath('/gpfs/ysm/home/lm978/meso_processing-master'));
%% input
%mainDir=fullfile('/gpfs/ysm/scratch60/cardin/lm978',Session); 
mainDir=fullfile('Z:\Imaging\Hannah\Raw data\MecP2',Session); 
%TmpFile='C:\Users\AHM\Documents\tmp.tif'; % temporary file in local ssd folder while loading images for speed, you can put this in C: but needs administrative access
TmpFile='temp'; 
outputFolder='X:\Hadas\Meso-imaging\CRISPR\traces_data\';
mkdir(strcat(outputFolder,'/',Session_after));
disp(strcat(Session,'processing')); 
%% user-selected input parameters 
params.Loc=location; %recording location, 'B449' or 'F238'
params.fsspike2=5000;
params.fsimaging=10; 
params.pupilSR=10; 
params.batchSize = 3000;
params.batches2load=10000;
params.regressUV=1; %option of regressing every pixel by UV 
params.moveLocal=1;%if you want to move the tiff file to a local drive before reading 
params.drawPlots=0;%plot all event channels and imaging channel, for test purposes
params.makeMovies=0;%make tif df/f movies for a segement of the session 
params.ImageFormat='tif';%choose between cxd or tif 
% visusal stim 
params.visStimAn=0; % 1 if vis stim are presented, 0 if not 
if strcmp(sessid,'Vis')
    params.visStimAn=1;
end 
params.airpuffAn=0; % if airpuffs or electrical stim are presented, indicate 1; 
if strcmp(sessid,'Air')
    params.airpuffAn=1;
end 
params.visDur=2;
params.visITI=5; 
%grabs/rcamp vs grabs only
params.signalsExtraction.firstCaFrame = 1;
params.signalsExtraction.sigs = color;% 'blueuv' or 'RCaMP_AC'
params.signalsExtraction.blueuvRatio = 1;
params.resizeRatio = 0.5;
% detrending and spatial filtering params
params.deterend.filtLen = 1000;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';
params.spatiotemporalfilter=1; %set it to 1 if you want to do spatial filtering  otherwise 0 

%run params
params.minRunDuration=5; %in seconds
params.minSitDuration=10; %in seconds 
params.ITITime=5; %  dead time after events in seconds  

%analysis window for events (airpuff and visstim)
params.preEventWin=5;
params.postEventWin=5;
params.baselineWin=2; %use this period to calculate baseline during within trial normalization 
%% choose right parameters depending on the recording location 
if strcmp(params.Loc,'B449')
params.angle2rotate=-180;
% spike2 channels for bcmm 
channels.BLUE = 1;%blue led 
channels.UV = 2;%uv led 
channels.FRAMETICKS = 8;%green/uv mesocam triggers
channels.PHOTO_DIODE = 4;%visual stim trigger
channels.WHEEL = 5;
channels.AIR_PUFF = 6; %this is either airpuff channel or electrical stim channel 
channels.PUPIL_CAMERA = 7;
channels.RED_MESO = 3;%red mesocam trigger
channels.GREEN=9;%green LED
channels.EEG = 10;%eeg continous signal 

elseif strcmp(params.Loc,'F238')
params.angle2rotate=90;
% spike2 channels for cardin 238
channels.BLUE = 1;%blue led 
channels.UV = 2;%uv led 
channels.FRAMETICKS = 3;%green/uv mesocam triggers
channels.PHOTO_DIODE = 4;%visual stim trigger
channels.WHEEL = 5;
channels.AIR_PUFF = 6; %this is either airpuff channel or electrical stim channel 
channels.PUPIL_CAMERA = 7;
channels.RED_MESO = 3;%red mesocam trigger
channels.EEG = 8;%eeg continous signal 
end 

if strcmp(cohort,'Cohort 4')
   params.angle2rotate=180;
end
%% get tiffs, smrx , visual csv path and run the main function 
tiffsPath = strcat(mainDir,'/meso');
smrxfilelist = (dir(fullfile(mainDir, '*.smrx')));
dataSmrxFile=fullfile(mainDir,smrxfilelist.name); 
dataVisStimList=(dir(fullfile(mainDir, '*.csv')));
dataVisStimFile=fullfile(mainDir,dataVisStimList.name); 
outputPath=fullfile(outputFolder,Session_after);

%% preprocessing step
disp('Starting pre-processing');
tic
main_generic_meso_processing_newv_C4_demix_align(tiffsPath, dataSmrxFile,dataVisStimFile,outputPath,TmpFile,channels,params);
toc
%% analysis step 1
%commented from tic to toc while we haven't had spike2 processed yet
%tic
% load preprocessed df/f if it doesnt exist 
%if ~exist('dFoF','var'),load(fullfile(outputPath, 'final_dFoF.mat'), 'dFoF', 'R', 'C');end 
%if ~exist('dFoF_parcells','var'),load(fullfile(outputPath, 'final_dFoF_parcels.mat'), 'dFoF_parcells');end 

%main_Level2_meso_processing(outputPath,params,dFoF,dFoF_parcells,R,C,tiffsPath) 
%toc 
end 


