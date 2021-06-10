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

function [dFoF,dFoF_parcells,R,C]= main_generic_meso_processing_newv_svdtest_14(tiffsPath, dataSmrxFile,dataVisStimFile, outputPath, TmpFile,channels,params)
addpath(genpath('../../../utils'));

cedpath = '../pre_processing_scripts/utils/CEDS64ML';
addpath(genpath('../parcellation/'));
addpath(genpath('../pre_processing_scripts/'));
addpath(genpath('..\regression'));
addpath(genpath('../correlation/'));
addpath(genpath('../meta_data_processing/'));

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
   visStimDur=params.visDur;
   visStimITI=params.visITI;
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
    % for Rcamp AC setting
    params.signalsExtraction.sigs = 'RCaMP_AC';
    params.signalsExtraction.blueuvRatio = 1;
    params.signalsExtraction.firstCaFrame = 1;
    resizeRatio = 0.5;
    
    % detrending params
    params.deterend.filtLen = 1000;
    params.deterend.filtcutoff = 0.001;
    params.deterend.method = 'FIR';
    angle2rotate = 180;
    spatiotemporalfilter=0; 
end

if exist(fullfile(tiffsPath, 'raw_mov.mat'), 'file')%exist(fullfile(outputPath, 'raw_mov.mat'), 'file')
   %load(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C');  
   load(fullfile(tiffsPath, 'raw_mov.mat'),'sigsMov','R','C');  %modified by lav
else
    %% load tiffs
    mkdir(outputPath);
    disp(['load tiffs from ' tiffsPath ]);
    if strcmp(imageFormat,'tif')
    mov = read_batches_tiffs(batchSize, tiffsPath, batches2load, resizeRatio,TmpFile,moveLocal);
    elseif strcmp (imageFormat,'cxd')
    mov=readCxdSingleFile(tiffsPath,resizeRatio);
    else
    error('Unknown file format'); 
    end 
    [R, C, ~] = size(mov);
    %% separating channels
    disp('Extract signals from tiffs');
    sigsMov = extract_signals_from_mov(mov, params.signalsExtraction); clear mov; 
    save(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C', '-v7.3');
end 

if exist(fullfile(outputPath, 'final_dFoF.mat'), 'file')
    load(fullfile(outputPath, 'final_dFoF.mat'), 'dFoF', 'R', 'C');   
    load(fullfile(outputPath, 'final_dFoF_parcels.mat'), 'dFoF_parcells');
else
    %% Registration
     C=size(sigsMov.blue,1)./R; 
    %if the data contains red-green-uv alternating frames
    if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
        disp(['Registering green and blue images for' tiffsPath]);
        %register blue and green images
        tformfile = fullfile(outputPath, 'tform_bluegreen.mat' );
        if ~exist(tformfile, 'file')
            [tform] = registerGreenBlue(sigsMov.green,sigsMov.blue,R,C,outputPath);%register two color images from two cameras
            save(tformfile, 'tform','R','C');
        else
            load(tformfile, 'tform','R','C');
        end
        sigsMov.green = transform_frames(sigsMov.green, tform, R, C);
        sigsMov.green=single(sigsMov.green); 
    end
     %% De-trend followed by df/f 
    names = fieldnames(sigsMov);
    ind1 = strcmp(names,'skippedframes' );
    ind2 = strcmp(names,'skippedchannels' );
   
    siginds = find(~ind1 & ~ind2); 
    
    for i = 1:length(siginds)
        sigsMov.(names{siginds(i)}) = imrotate_vectored(sigsMov.(names{siginds(i)}),R, C, params.angle2rotate);%rotate the movie so the front of the brain is at the bottom 
        firstFrame.(names{siginds(i)})=(sigsMov.(names{siginds(i)})(:,1000)); %extract one frame from the raw movie 
        firstFrame.(names{siginds(i)})=imadjust(mat2gray( firstFrame.(names{siginds(i)}))); %this first frame is for plotting the parcel overlay in the parcellation step, looks better than df/f for identifying anatomical landmarks        
        h=figure; subplot(1,2,1),imagesc(reshape(firstFrame.(names{siginds(i)}),R,C));title('Before'); %draw a figure before removal of saturated pixels        
        maxPixelValue=(max(sigsMov.(names{siginds(i)}),[],2));
        saturatedPixelIdx{i}=find(maxPixelValue>=65535); %find saturated pixels 
        sigsMov.(names{siginds(i)})(saturatedPixelIdx{i},:)=0;  %remove any saturated pixel 
        subplot(1,2,2),imagesc(reshape(sigsMov.(names{siginds(i)})(:,1000),R,C));%draw a figure after removal of saturated pixels
        runTitle=strcat('After',' number of saturated Pixels = ',num2str(numel(saturatedPixelIdx{i})));
        title(runTitle); imageName=fullfile(outputPath,strcat('RemovalSaturatedPixels',names{i})); saveas(h,imageName);        
        disp(['Detrending ' names{siginds(i)}]);
        [sigsMov.(names{siginds(i)}),sigsBaseline.(names{siginds(i)})]=detrend_all_pixels(sigsMov.(names{siginds(i)}), params.deterend);
        disp(['Extracting dFoF ' names{siginds(i)}]);
        dFoF.(names{siginds(i)}) = sigsMov.(names{siginds(i)})./ sigsBaseline.(names{siginds(i)});%calculated df/f with f0 as the low pass filtered signal 
       
    end
    clear mov sigsMov maxPixelValue     
    %% register movies to allen atlas and parcellate into brain regions 
    disp(['Atlas Registration for' tiffsPath]);        
    load('parcells_updated121519.mat'); parcells=parcells_new;parcells_template=(mat2gray(parcells.CombinedParcells));%load new parcells    
    tformfile = fullfile(tiffsPath, 'tform_blue.mat' );
    if ~exist(tformfile, 'file')
         [tform,R,C] = get_alignment_CotrolPoints_transform(firstFrame.blue,parcells_template,R,C,parcells.movingPoints,parcells.fixedPoints); % do alignment based on control points manually selected on GUI, you can pass what template you want to use
        save(tformfile, 'tform','R','C');
    else
        load(tformfile, 'tform','R','C');
    end
    
    names = fieldnames(dFoF);
    for i = 1:length(names)
        dFoF.(names{i}) = transform_frames(dFoF.(names{i}), tform, R, C);
        firstFrame_t.((names{i}))=transform_frames(firstFrame.(names{i}), tform, R, C);
        [h1]=plot_parcell_overlay(firstFrame_t.(names{i}),R,C,1,parcells.indicators);
        if ~isempty(h1)
            imageName1=fullfile(outputPath,strcat('Parcels_',names{i}));
            saveas(h1,imageName1);
        end
    end
   
    %% Spatial SVD method for each color pair 
    if params.regressUV  %do the SVD method if set to 1 otherwise you skip this UV artifact removal 
    %parpool(feature('NumCores'))
    % parpool('threads')
    % load brain mask
    load('brainMask.mat','mask');
    maskinds = find(mask);
    naninds=find(isnan(dFoF.blue(:,1)));
    for i=1:length(naninds)
        [naninds_sub(i,1),naninds_sub(i,2)]=ind2sub([R,C],naninds(i));
        mask(naninds_sub(i,1),naninds_sub(i,2))=0;
    end
    %ensure that all fields have the same number of elements 
    names = fieldnames(dFoF);
    minsize=size(dFoF.uv,2); 
    for i =1:length(names)
         minsize=min(size(dFoF.(names{i}),2),minsize); 
    end 
    %smooth uv first 
    dFoF.uv=smoothdata(dFoF.uv,2,'gaussian',5); 
    dFoF.uv=dFoF.uv(:,1:minsize);
    %save(fullfile(outputPath, 'Blue_UV_BeforeSVD'),'dFoF.uv', 'dFoF.blue' );
    %run spatial svd with a patch size of 11
    for i = 1:length(names)
        if  ~strcmp(names{i},'uv')
            if exist(fullfile(outputPath, strcat('SVDOut',names{i},'.mat')), 'file') % if the output already exists load that
                disp(strcat('LoadingSVDOutput ',names{i}));
                load(fullfile(outputPath, strcat('SVDOut',names{i})),'output_sig');
                dFoF.(names{i})=output_sig;clear output_sig;
            else
                disp(strcat('SVD correction processing',names{i}));
                dFoF.(names{i})=dFoF.(names{i})(:,1:minsize);
                [output_sig, sig_blue, sig_uv] = spatial_regression_fast(dFoF.(names{i}), dFoF.uv, mask, 11, 0);
                save(fullfile(outputPath, strcat('SVDOut',names{i})),'output_sig', 'sig_blue','sig_uv','-v7.3' );
                dFoF.(names{i})=output_sig;clear output_sig;
            end
        end
    end
    end
  %% parcellate data   
  for i = 1:length(names)  
       disp(['Parcellating' names{i}]); 
       dFoF_parcells.(names{i}) = pixels2rois(dFoF.(names{i}), parcells); 
  end 
   
  %% make a movie for a segment of df/f 
  if params.makeMovies
      timeSeg=1:2000;
      for i = 1:length(names)
          currMov=mat2gray(reshape(dFoF.(names{i})(:,timeSeg),R,C,length(timeSeg)));
          outFile=fullfile(outputPath,strcat('Dff_Movie_',names{i},num2str(timeSeg(1)),'-',num2str(timeSeg(end)),'.tif'));
          h1=currMov(:,:,timeSeg(1));
          imwrite(im2uint16(h1),outFile ,'tif');
          for k=timeSeg(2):timeSeg(end)
              h1 = currMov(:,:,k);
              imwrite(im2uint16(h1), outFile,'WriteMode', 'append');
          end
      end
  end
  save(fullfile(outputPath, 'final_dFoF.mat'), 'dFoF','R','C','firstFrame', '-v7.3');
  save(fullfile(outputPath, 'final_dFoF_parcels.mat'), 'dFoF_parcells', '-v7.3');  
end 

%% extract timestamps from smrx file and align imaging data to spike2 data
disp(['Extractign smrx timestamps for' tiffsPath]);
if exist(fullfile(outputPath, 'smrx_signals.mat'), 'file')
    if params.visStimAn
        load(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx','params');
    else
        load(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps');
    end
else
    display(strcat('loading spike 2 file: ',dataSmrxFile));
    [timing, channels_data] = process_spike2_two_cams(cedpath, outputPath,dataSmrxFile, fsspike2,channels,minRunDuration,minSitDuration,ITITime,visStimDur,visStimITI,pupilSR);
    
    % get timestamps for all relevant events
    timing.mesostart=timing.mesostart(1:length(timing.mesoend)); % remove the last mesostart if there is an extra one without a corresponding end
    timestamps.timaging = (timing.mesostart+timing.mesoend)/2;
    
    switch params.signalsExtraction.sigs
        case 'blueuv'
            timestamps.timaging=timestamps.timaging(1:2:end);
        case 'RCaMP_AC'
            timestamps.timaging=timestamps.timaging(1:3:end);
    end
    
    timestamps.timaging=timestamps.timaging(params.deterend.filtLen/2:end);%remove the first filtLen/2 timestamps because there is a filtering artifact and those samples are removed in the detrending step  
    timestamps.timaging=timestamps.timaging(1:size(dFoF_parcells.blue,2));%remove excess timestamps if you don't have corresponding frames
    timestamps.visStim=timing.stimstart; 
    timestamps.airpuff=timing.airpuffstart(timing.airpuffstart>(timestamps.timaging(1)+ preEventWin) & timing.airpuffstart<(timestamps.timaging(end)-postEventWin));%remove events if they occur outside imaging period 
    timestamps.wheelOn=timing.wheelOn(timing.wheelOn>(timestamps.timaging(1)+ minSitDuration) & timing.wheelOn<(timestamps.timaging(end)-minRunDuration));%remove events if they occur outside imaging period 
    
    %% if you want to analyze visual stim
    if params.visStimAn
        visData=xlsread(dataVisStimFile);
        contrasts=visData(:,1);
       
        if numel(contrasts)==numel(timestamps.visStim)%get timestamps for each contrast
           visStimidx= find(timing.stimstart>(timestamps.timaging(1)+preEventWin) & timing.stimstart<(timestamps.timaging(end)-postEventWin));
           timestamps.visStim =timing.stimstart(visStimidx);%remove events if they occur outside of imaging period
           contrasts=contrasts(visStimidx); %remove events if they occur outside of imaging period
           uniqContrasts=unique(contrasts);
            for i=1:length(uniqContrasts)
                contrast_Idx{uniqContrasts(i)}=find(contrasts==uniqContrasts(i));
                timestamps.contrasts{uniqContrasts(i)}=timestamps.visStim(contrast_Idx{uniqContrasts(i)});
            end
        else
            error('error:number of stims in excel and spike2 do not mtatch')
        end
        
        save(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','uniqContrasts','contrast_Idx','params');
    else
        save(fullfile(outputPath, 'smrx_signals.mat'), 'timing', 'channels_data','timestamps','params');
    end
end
%% draw ROI around V1 and plot all event channels 
if params.drawPlots

    switch params.signalsExtraction.sigs
        case 'blueuv'
            numPlots=6;
            meanFrame=zeros(R,C);
            if params.visStimAn
            for k = 3:length(timestamps.visStim)
                ind_s = findClosestDouble(timestamps.timaging, timestamps.visStim(k));
                ind_e = findClosestDouble(timestamps.timaging, timestamps.visStim(k)+2);
                meanFrame = meanFrame + reshape(mean(dFoF.blue(:, ind_s:ind_e),2), R, C);
            end
            elseif params.airpuffAn
            for k = 3:length(timestamps.airpuff)
                ind_s = findClosestDouble(timestamps.timaging, timestamps.airpuff(k));
                ind_e = findClosestDouble(timestamps.timaging, timestamps.airpuff(k)+2);
                meanFrame = meanFrame + reshape(mean(dFoF.blue(:, ind_s:ind_e),2), R, C);
            end
            end 
            t_spike2 = linspace(0, (length(channels_data.mesoframe)-1)/fsspike2, length(channels_data.mesoframe));
            figure;imagesc(meanFrame)
            title('Please select V1');
            drawmask= drawpolygon(gca);
            mask=drawmask.createMask;
            imagesc(mask);colorbar;
            v1timetrace_Blue = mean(dFoF.blue(find(mask==1), :));
            v1timetrace_UV = mean(dFoF.uv(find(mask==1), :));
        case 'RCaMP_AC'
            numPlots=7;
            meanFrame=zeros(R,C);
            if params.visStimAn
            for k = 3:length(timestamps.visStim)
                ind_s = findClosestDouble(timestamps.timaging, timestamps.visStim(k));
                ind_e = findClosestDouble(timestamps.timaging, timestamps.visStim(k)+2);
                meanFrame = meanFrame + reshape(mean(dFoF.green(:, ind_s:ind_e),2), R, C);
            end
            elseif params.airpuffAn
            for k = 3:length(timestamps.airpuff)
                ind_s = findClosestDouble(timestamps.timaging, timestamps.airpuff(k));
                ind_e = findClosestDouble(timestamps.timaging, timestamps.airpuff(k)+2);
                meanFrame = meanFrame + reshape(mean(dFoF.green(:, ind_s:ind_e),2), R, C);
            end
            end 
            t_spike2 = linspace(0, (length(channels_data.mesoframe)-1)/fsspike2, length(channels_data.mesoframe));
            figure;imagesc(meanFrame)
            title('Please select V1');
            drawmask= drawpolygon(gca);
            mask=drawmask.createMask;
            imagesc(mask);colorbar;
            v1timetrace_Green = mean(dFoF.green(find(mask==1), :));
            v1timetrace_Blue = mean(dFoF.blue(find(mask==1), :));
            v1timetrace_UV = mean(dFoF.uv(find(mask==1), :));
    end
    figure;
    h(1) = subplot(numPlots,1,1);
    plot(t_spike2,channels_data.mesoframe);
    title('mesoframes');
    h(2) = subplot(numPlots,1,2);
    plot(t_spike2,channels_data.diode)
    title('stim');
    h(3) = subplot(numPlots,1,3);
    plot(t_spike2,channels_data.air_puff)
    title('air puff');
    h(4) = subplot(numPlots,1,4);
    plot(t_spike2, channels_data.wheelspeed);
    title('wheelspeed');
    h(5) = subplot(numPlots,1,5);
    plot(timestamps.timaging, v1timetrace_Blue);
    ylim([-0.02 0.05]);
    title('Blue');
    h(6) = subplot(numPlots,1,6);
    plot(timestamps.timaging, v1timetrace_UV);
    ylim([-0.02 0.05]);
    title('UV');
    if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
    h(7) = subplot(numPlots,1,7);
    plot(timestamps.timaging, v1timetrace_Green);
    title('Green');
    ylim([-0.02 0.05]);
    end 
    xlabel('Time [secs]');
    linkaxes(h,'x');
    xlim([timestamps.timaging(1) timestamps.timaging(end)]);

%% dff for selected ROI only around visual stimulus 
if params.visStimAn
h=figure; 
for k=1:length(uniqContrasts)
currContrast=uniqContrasts(k); 
eventTS=timestamps.contrasts{1,currContrast};preEventWin=2;postEventWin=5; eventTS=eventTS(:)';
[trial_dff,timeStamps]= TrialTimeArrangeDff(v1timetrace_Green,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
%trial_dff=fillmissing(trial_dff,'linear',1,'EndValues','nearest');%interpolate
plot(mean(timeStamps,2),nanmean(trial_dff,2)); hold on; ylabel('Dff'),xlabel('Time(s)'); 
end 
legend('2%','5%','7%','20%','50%','100%'); 
saveas(h,fullfile(outputPath,'DrawnROIDFF_green'));
end 
end 
%get dff from area in left primary visual cortex wtih maximal activation
%during 100% contrast presentation 
if params.visStimAn && ~exist(fullfile(outputPath, 'maxVisualTrace.mat'), 'file')
if ~exist('parcells','var'), load('parcells_updated121519.mat'); parcells=parcells_new; end 
V1parcell=parcells.indicators (:,:,2); %mask for left visual cortex 
V1parcellIdx=find(V1parcell==1); 
visActivity_time=(dFoF.green(V1parcellIdx, :));
eventTS=timestamps.contrasts{1,100};preEventWin=0;postEventWin=2; eventTS=eventTS(:)';
[trial_dff,timeStamps]= TrialTimeArrangeDff(visActivity_time,timestamps.timaging,fsimaging,preEventWin,eventTS,postEventWin);
visActivity_trial=mean(trial_dff,1);
visActivity=squeeze(mean(visActivity_trial,2));
[Ms,Sortidx] = sort(visActivity,'descend');% Sort asending along time dimension
sortedParcellIdx=V1parcellIdx(Sortidx); 
maxVisSortedIdx=sortedParcellIdx(1:ceil(length(Ms)*0.2));%top 20% of values 
h1=figure; tmp=mean(dFoF.green(:, 1:1000),2);imageV=reshape(tmp, R, C); % plot the top 20%
subplot(2,1,1);imagesc(imageV);tmp(maxVisSortedIdx)=1000; top20ImageV=reshape(tmp,R,C);  hold on;
subplot(2,1,2);imagesc(top20ImageV); hold off;
saveas(h1,fullfile(outputPath,'top20Pixels_visStim'));
save(fullfile(outputPath, 'maxVisualTrace.mat'),'maxVisSortedIdx');  
end 
end


