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

function main_generic_meso_processing_newv_C4_demix_align(tiffsPath, dataSmrxFile,dataVisStimFile, outputPath, TmpFile,channels,params)
addpath(genpath('../../../utils'));

cedpath = 'X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\pre_processing_scripts\utils\CEDS64ML';
addpath(genpath('..\s_meso_processing-master'));
addpath(genpath('..\meta_data_processing'));
addpath(genpath('../correlation/'));
addpath(genpath('../s_meso_processing-master/pre_processing_scripts/'));

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

if exist(fullfile(outputPath, 'raw_mov.mat'), 'file')%exist(fullfile(outputPath, 'raw_mov.mat'), 'file')
   %load(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C');  
   load(fullfile(outputPath, 'raw_mov.mat'),'sigsMov','R','C');  %modified by lav
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
    %PixxTime_bl=sigsMov.blue;
    %PixxTime_uv=sigsMov.uv;
    %save(fullfile(outputPath, 'raw_pixxtime.mat'),'PixxTime_bl','PixxTime_uv','R','C', '-v7.3');
end 
if exist(fullfile(outputPath, 'final_dFoF.mat'), 'file')
    load(fullfile(outputPath, 'final_dFoF.mat'), 'dFoF', 'R', 'C');   
    load(fullfile(outputPath, 'final_dFoF_parcels.mat'), 'dFoF_parcells');
else
    %% Registration
    C=size(sigsMov.blue,1)./R;
    %if the data contains red-green-uv alternating frames
%     if strcmp(params.signalsExtraction.sigs,'RCaMP_AC')
%         disp(['Registering green and blue images for' tiffsPath]);
%         %register blue and green images
%         tformfile = fullfile(tiffsPath, 'tform_bluegreen.mat' );
%         if ~exist(tformfile, 'file')
%             [tform] = registerGreenBlue(sigsMov.green,sigsMov.blue,R,C,outputPath);%register two color images from two cameras
%             save(tformfile, 'tform','R','C');
%         else
%             load(tformfile, 'tform','R','C');
%         end
%         sigsMov.green = transform_frames(sigsMov.green, tform, R, C);
%         sigsMov.green=single(sigsMov.green);
%     end
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
%         maxPixelValue=(max(sigsMov.(names{siginds(i)}),[],2));
%         saturatedPixelIdx{i}=find(maxPixelValue>=65535); %find saturated pixels 
%         sigsMov.(names{siginds(i)})(saturatedPixelIdx{i},:)=0;  %remove any saturated pixel 
%         subplot(1,2,2),imagesc(reshape(sigsMov.(names{siginds(i)})(:,1000),R,C));%draw a figure after removal of saturated pixels
%         runTitle=strcat('After',' number of saturated Pixels = ',num2str(numel(saturatedPixelIdx{i})));
%         title(runTitle); imageName=fullfile(outputPath,strcat('RemovalSaturatedPixels',names{i})); saveas(h,imageName);        
%         disp(['Detrending ' names{siginds(i)}]);
%         [sigsMov.(names{siginds(i)}),sigsBaseline.(names{siginds(i)})]=detrend_all_pixels(sigsMov.(names{siginds(i)}), params.deterend);
%         disp(['Extracting dFoF ' names{siginds(i)}]);
%         %dFoF.(names{siginds(i)}) = sigsMov.(names{siginds(i)})./ sigsBaseline.(names{siginds(i)});%calculated df/f with f0 as the low pass filtered signal 
       
    end
    clear mov sigsMov maxPixelValue       

    disp(['Atlas Registration for' tiffsPath]);
    load('parcells_updated121519.mat'); parcells=parcells_new;parcells_template=(mat2gray(parcells.CombinedParcells));%load new parcells
    tformfile = fullfile(outputPath, 'tform_blue.mat' );
    if ~exist(tformfile, 'file')
        [tform,R,C] = get_alignment_CotrolPoints_transform(firstFrame.blue,parcells_template,R,C,parcells.movingPoints,parcells.fixedPoints); % do alignment based on control points manually selected on GUI, you can pass what template you want to use
        save(tformfile, 'tform','R','C');
    else
    end
    names = {'blue','uv'};
    for i = 1:length(names)
        firstFrame_t.((names{i}))=transform_frames(firstFrame.(names{i}), tform, R, C);
        [h1]=plot_parcell_overlay(firstFrame_t.(names{i}),R,C,1,parcells.indicators);
        if ~isempty(h1)
            imageName1=fullfile(outputPath,strcat('Parcels_',names{i}));
            saveas(h1,imageName1);
        end
    end
    disp('done')
end
end


