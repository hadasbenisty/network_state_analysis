clear;close all;
restoredefaultpath;
addpath(genpath('../network_state_analysis'))
addpath(genpath('../network_state_analysis/meta_data_processing'))

%addpath('X:\Hadas\Meso-imaging\Antara\results\code\functions');

csvfile = '../network_state_analysis/meta_data_processing/antara_exp/Processing_Pipeline_Full_expression_andHannah.xlsx';
processing_pipeline=readtable(csvfile);

%preprocess
processing_pipeline=processing_pipeline{strcmp(processing_pipeline{:,6},'Cohort 4'),1:7}
airindex=zeros(size(processing_pipeline,1),1);sponindex=zeros(size(processing_pipeline,1),1);

for i=1:size(processing_pipeline,1)
    if strcmp(char(processing_pipeline{i,4}),'Air')
        airindex(i)=1;
    elseif strcmp(char(processing_pipeline{i,4}),'Spon')
        sponindex(i)=1;
    else
    end
end

%% pupil manual step
for j=2
    if j==1
        index=airindex;sessions_it=processing_pipeline(find(airindex==1),2);sessid='Air';
    elseif j==2
        index=sponindex;sessions_it=processing_pipeline(find(sponindex==1),2);sessid='Spon';
    end
    for i=1:size(find(index==1),1)
        if contains(string(processing_pipeline{i,1}),'Emx')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Emx_');
        elseif contains(string(processing_pipeline{i,1}),'Control')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Control_');
        elseif contains(string(processing_pipeline{i,1}),'Lhx6')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Lhx6_');
        end
        Session=strcat(string(processing_pipeline{i,7}),'_',string(processing_pipeline{i,1}));
        if strcmp(Session,'emx')
            Session=strrep(Session,'emx','Emx');
        elseif strcmp(Session,'ctrl')
            Session=strrep(Session,'ctrl','Control');
        elseif strcmp(Session,'lhx6')
            Session=strrep(Session,'lhx6','Lhx6');
        end

        tiffsPath=char(fullfile('Z:\Imaging\Hannah\Raw data\MecP2\',string(processing_pipeline{i,1})));
        outputFolder='X:\Hadas\Meso-imaging\CRISPR\traces_data';
        disp(strcat('Processing ',Session));
        
        [params,channels] = get_channels_param(sessid,'F238','blueuv');
        %% get tiffs, smrx , visual csv path and run the main function
        smrxfilelist = (dir(fullfile(tiffsPath, '*.smrx')));
        dataSmrxFile=fullfile(tiffsPath,smrxfilelist.name);
        outputPath=fullfile(outputFolder,Session);
        mkdir(char(strcat(outputPath)));
        %% preprocessing step
        tic
         %if (i~=12)&(i~=3)
        %images2Avi(strcat(tiffsPath,'\pupil'), 'pupil', 9)
        detect_axis_for_pupil(outputPath, strcat(tiffsPath,'\pupil'), false)
         %else
         %end
        toc
        clearvars -except i j sessions_it airindex index sponindex sessid processing_pipeline
    end
end


%% pupil automatic step
for j=2
    if j==1
        index=airindex;sessions_it=processing_pipeline(find(airindex==1),2);sessid='Air';
    elseif j==2
        index=sponindex;sessions_it=processing_pipeline(find(sponindex==1),2);sessid='Spon';
    end
    for i=1:size(find(index==1),1)
        if contains(string(processing_pipeline{i,1}),'Emx')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Emx_');
        elseif contains(string(processing_pipeline{i,1}),'Control')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Control_');
        elseif contains(string(processing_pipeline{i,1}),'Lhx6')
            processing_pipeline{i,1}=erase(string(processing_pipeline{i,1}),'Lhx6_');
        end
        Session=strcat(string(processing_pipeline{i,7}),'_',string(processing_pipeline{i,1}));
        if strcmp(Session,'emx')
            Session=strrep(Session,'emx','Emx');
        elseif strcmp(Session,'ctrl')
            Session=strrep(Session,'ctrl','Control');
        elseif strcmp(Session,'lhx6')
            Session=strrep(Session,'lhx6','Lhx6');
        end
        tiffsPath=char(fullfile('Z:\Imaging\Hannah\Raw data\MecP2\',string(processing_pipeline{i,1})));
        outputFolder='X:\Hadas\Meso-imaging\CRISPR\traces_data';
        disp(strcat('Processing ',Session));
        
        [params,channels] = get_channels_param(sessid,'F238','blueuv');
        %% get tiffs, smrx , visual csv path and run the main function
        tiffsPath = tiffsPath;
        smrxfilelist = (dir(fullfile(tiffsPath, '*.smrx')));
        dataSmrxFile=fullfile(tiffsPath,smrxfilelist.name);
        outputPath=fullfile(outputFolder,Session);
        mkdir(char(strcat(outputPath)));
        %% preprocessing step
        tic
        %step 3
        % eyeblinkdetect-semi automatic
        % semi_auto_eye_blink(outputpath, daysession, animal, ANALOG_RATE, overwrite_binary2_blink_blinksummary);
        
        % spontaneous eyeblink detect, binarize whisking and normalize pupil size
        % detects the timing of spont blink
        %detect_timing_spont_blink(animal, daysession, outputpath, overwirte_binary_roi_2_spontblink);
        if exist(fullfile(outputPath,'\roi.mat'), 'file')==0 &&exist(fullfile(outputPath,'\pupil.mat'), 'file')==0
            % normalize pupil size based on dark-light changes/minimum pupil/eye size
            try
                extract_roi_intensity_detect_pupil_size_and_whiskering(strcat(tiffsPath,'\pupil'), outputPath, false)
                normalize_pupil_by_dark_light_changes(outputPath, tiffsPath,false)
            catch
                disp(strcat('ERROR',outputPath))
            end
        end
        toc
        clearvars -except i j sessions_it airindex index sponindex sessid processing_pipeline
    end
end

% for j=1
%     if j==1
%         index=airindex;sessions_it=processing_pipeline(find(airindex==1),2);sessid='Air';
%     elseif j==2
%         index=sponindex;sessions_it=processing_pipeline(find(sponindex==1),2);sessid='Spon';
%     end
%     for i=1:size(find(index==1),1)
%         Session=strcat(string(processing_pipeline{i,5}),'_',string(processing_pipeline{i,1}));
%         mainDir=char(fullfile('X:\CardinLab\Antara\',Session));
%         outputFolder='X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData';
%         disp(strcat('Processing ',Session));
%         [params,channels] = get_channels_param(sessid,'F238','blueuv');
%         %% get tiffs, smrx , visual csv path and run the main function
%         tiffsPath = mainDir;
%         smrxfilelist = (dir(fullfile(mainDir, '*.smrx')));
%         dataSmrxFile=fullfile(mainDir,smrxfilelist.name);
%         outputPath=fullfile(outputFolder,Session);
%         mkdir(char(strcat(outputPath)));
%         %% preprocessing step
%         tic
%         disp(strcat('processing',outputPath))
%         main_generic_meso_processing_newv_svdtest_spike2(tiffsPath, dataSmrxFile,'temp',outputPath,'temp',channels,params);
%         toc
%         clearvars -except i j sessions_it airindex index sponindex sessid processing_pipeline
%     end
% end