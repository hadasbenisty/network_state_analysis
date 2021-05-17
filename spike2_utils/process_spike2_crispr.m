

function process_spike2_crispr(cedpath, tiffsPath, dataSmrxFile,outputPath, channels,params,Session)

%smrx_signals_v3 is with 0.5sec window smoothing wheel
%smrx_signals_v4 is with 1sec window smoothing wheel
if isfile(fullfile(outputPath, 'smrx_signals_v4.mat'))
    disp('file exist');
    return;
end

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

% detrending params
params.deterend.filtLen = 150;
params.deterend.filtcutoff = 0.001;
params.deterend.method = 'FIR';

%% extract timestamps from smrx file and align imaging data to spike2 data
disp(['Extractign smrx timestamps for' tiffsPath]);

display(strcat('loading spike 2 file: ',dataSmrxFile));
fsspike2=5000;

[timing, channels_data] = process_spike2_two_cams(cedpath, outputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,pupilSR,tiffsPath);

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


mkNewDir(outputPath)
save(char(strcat(outputPath, '\smrx_signals_v4.mat')), 'timing', 'channels_data','timestamps','params');
end



