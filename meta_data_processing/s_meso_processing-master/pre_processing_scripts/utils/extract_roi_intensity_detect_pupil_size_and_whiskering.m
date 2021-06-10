function extract_roi_intensity_detect_pupil_size_and_whiskering(inputpath, outputpath, overwrite)
% Extract roi intensity detect pupil size_and whiskering.
% Reads mp4/avi file and _axisroi.mat file creates output files eding
%  with _roi, _pupil, _whisker
%
% input -
%   datapath                       - path to read input files
%   outputpath                     - path to save output files
%   daysession                     - a list of days for processing
%   animal                         - a list of animals as a char array
%   overwrite_avimp4__axisroi_2_roi_pupil_whisker - if to overwrite output
%
%   Written by Lan Tang
%   editted by Hadas Ben Esti 6/28/19
addpath(genpath('../pre_processing_scripts/utils'));
pointdist=@(a,b) sqrt(sum((a-b).^2));
warning off;
data_time_stamp_filename=fullfile(inputpath);
output_time_stamp_filename=fullfile(outputpath);
% 
% if contains(output_time_stamp_filename,'Control')
%     Session = strrep(output_time_stamp_filename,'9 C','9_C');
% elseif contains(output_time_stamp_filename,'Lhx6')
%     Session = strrep(pre_Session,'6 C','6_C');
% elseif contains(output_time_stamp_filename,'Emx')
%     Session = strrep(output_time_stamp_filename,'x C','x_C');
% end
if exist(strcat(output_time_stamp_filename,'_roi.mat'), 'file') &&...
        exist(strcat(output_time_stamp_filename,'_pupil.mat'), 'file') &&...
        exist(strcat(output_time_stamp_filename,'_whisker.mat'), 'file') && ...
        ~overwrite
else
    display(strcat('pupil and blink detecting: ',data_time_stamp_filename));
    video=dir(fullfile(data_time_stamp_filename,'*.mp4'));
    if isempty(video)
        video=dir(fullfile(data_time_stamp_filename,'*.avi'));
    end
    mov=fullfile(video.folder,video.name);
    
    if exist(fullfile(output_time_stamp_filename,'axisroi.mat'),'file')
        load(fullfile(output_time_stamp_filename,'axisroi.mat'))
        
        %[roiint,areaii,centerx,centery,framenum]=lt_circlesdetect(mov,longaxis,roi,0.995);
        %change this
        [roiint,areaii,centerx,centery,circularity,framenum]=pupildetect_bw(mov,longaxis,roi,imagthres,pupilthres);
        
        %[roiint,areaii,centerx,centery,circularity,whisker,framenum]=pupwhiskdetect_mesoimaging(mov,eyerect,facerect,roi_eye,roi_face,thres,thres2,thres3);
        
        %roiint=sgolayfilt(roiint,1,91);
        save(fullfile(output_time_stamp_filename,'\roi.mat'),'roiint');
        save(fullfile(output_time_stamp_filename,'\pupil.mat'),'areaii','centerx','centery','circularity');
        clearvars -except filenum daysession animal pathname;
    else
        display(strcat(data_time_stamp_filename,'_axisroi.mat not exist'));
    end
end
end
