%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main script for processing meso imaging
%
% Input
%   - tiffsPath      - path to tiffs
%   - data_smr_file  - path to smrx file
%   - outputPath     - output path
%   - params         - a struct with fields:
%                       * fsspike2 - sampling freq of spike2 (default =
%                       5e3Hz)
%                       * fsimaging - sampling freq of imaging (default =
%                       30Hz)
% Written by Hadas Benisty
% 11/15/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_generic_meso_processing_newv_svdtest_spike2(tiffsPath, dataSmrxFile,dataVisStimFile, outputPath, TmpFile,channels,params)
addpath(genpath('../../../utils'));

cedpath = 'X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\pre_processing_scripts\utils\CEDS64ML';
addpath(genpath('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master'));
% addpath(genpath('../pre_processing_scripts/'));
% addpath(genpath('..\regression'));
% addpath(genpath('../correlation/'));
% addpath(genpath('../meta_data_processing/'));

if nargin == 7
    fsspike2 = params.fsspike2;
    fsimaging = params.fsimaging;
    pupilSR=params.pupilSR; 
    batchSize = params.batchSize;  
    resizeRatio=params.resizeRatio;
    angle2rotate=params.angle2rotate; 
    batches2load=params.batches2load;
    spatiotemporalfilter=params.spatiotemporalfilter;
    moveLocal=params.moveLocal; 
   imageFormat=params.ImageFormat;
   regressUV= params.regressUV;
   minRunDuration=params.minRunDuration;
   minSitDuration=params.minSitDuration;
   ITITime=params.ITITime; 
   preEventWin=params.preEventWin;
   postEventWin=params.postEventWin;
   baselineWin=params.baselineWin; 
else
    fsspike2 = 5e3;
    fsimaging=10; 
    batchSize = 3000;
    batches2load=1000;
    moveLocal=1; 
    imageFormat='tif';
    regressUV=0; 
    % for blue-uv settings
    params.signalsExtraction.firstCaFrame = 1;
    params.signalsExtraction.sigs = 'blueuv';
    params.signalsExtraction.blueuvRatio = 1;
    resizeRatio = 1;
    
    % detrending params
    params.deterend.filtLen = 150;
    params.deterend.filtcutoff = 0.001;
    params.deterend.method = 'FIR';
    angle2rotate = 180;
    spatiotemporalfilter=0; 
end
%% extract timestamps from smrx file and align imaging data to spike2 data
disp(['Extractign smrx timestamps for' tiffsPath]);
if 0%exist(fullfile(char(outputPath), 'smrx_signals_v2.mat'), 'file')
    if params.visStimAn
        load(fullfile(outputPath, 'smrx_signals_v2.mat'), 'timing', 'channels_data','timestamps');
    else
        load(fullfile(outputPath, 'smrx_signals_v2.mat'), 'timing', 'channels_data','timestamps');
    end
else
    
    display(strcat('loading spike 2 file: ',dataSmrxFile));
    fsspike2=5000;
    
    [timing, channels_data] = process_spike2_two_cams(cedpath, outputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,pupilSR,tiffsPath)
    
    timingwheelOn=timing.allwheelon;
    % get timestamps for all relevant events
    timing.mesostart=timing.mesostart(1:length(timing.mesoend)); % remove the last mesostart if there is an extra one without a corresponding end
    timestamps.timaging = (timing.mesostart+timing.mesoend)/2;
    
    switch params.signalsExtraction.sigs
        case 'blueuv'
            timestamps.timaging=timestamps.timaging(1:2:end);
        case 'RCaMP_AC'
            timestamps.timaging=timestamps.timaging(1:3:end);
    end
    
    %timestamps.timaging=timestamps.timaging(params.deterend.filtLen/2:end);%remove the first filtLen/2 timestamps because there is a filtering artifact and those samples are removed in the detrending step  
    %timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));%remove excess timestamps if you don't have corresponding frames
   % timestamps.visStim=timing.stimstart; 
 
    timestamps.airpuff=timing.airpuffstart(timing.airpuffstart>(timestamps.timaging(1)+ preEventWin) & timing.airpuffstart<(timestamps.timaging(end)-postEventWin));%remove events if they occur outside imaging period 
    timestamps.wheelOn=timing.wheelOn(timing.wheelOn>(timestamps.timaging(1)+ minSitDuration) & timing.wheelOn<(timestamps.timaging(end)-minRunDuration));%remove events if they occur outside imaging period 
    
    
if contains(outputPath,'Control')
    outputPath = strrep(outputPath,'9 C','9_C');
elseif contains(outputPath,'Lhx6')
    outputPath = strrep(outputPath,'6 C','6_C');
elseif contains(outputPath,'Emx')
    outputPath = strrep(outputPath,'x C','x_C');
end
    save(char(strcat(outputPath, '\smrx_signals_v2.mat')), 'timing', 'channels_data','timestamps','params');
end
%process pupil
end


