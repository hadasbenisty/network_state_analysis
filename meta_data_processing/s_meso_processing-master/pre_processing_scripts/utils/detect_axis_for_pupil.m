function detect_axis_for_pupil(outputpath, inputpath, overwrite)
% Semi-supervised for detection of pupil axis. Reads mp4 or avi file
%  and creates an output file with a '_axisroi' suffix
% input -
%   datapath                       - path to read input files
%   outputpath                     - path to save output files
%   daysession                     - a list of days for processing
%   animal                         - a list of animals as a char array
%   overwrite_avimp4_2axisroi      - if to overwrite output file if exists
%   Written by Lan Tang
%   editted by Hadas Ben Esti 6/28/19
addpath(genpath('../pre_processing_scripts/utils'));
outfile = fullfile(outputpath, 'axisroi.mat');
if exist(outfile,'file')&& ~overwrite
else
    data_time_stamp_filename=fullfile(inputpath);
    display(strcat('draw long axis and roi:  ', data_time_stamp_filename));
    video=dir(fullfile(data_time_stamp_filename,'*.mp4'));
    
    if isempty(video)
        video=dir(fullfile(data_time_stamp_filename,'*.avi'));
    end
    
    mov=fullfile(video.folder,video.name);

    if exist(mov,'file')
        close all;
        %[roi_eye,roi_face,eyerect,facerect,thres,thres2,thres3]=longaxisdraw_meso_new(mov);
        [longaxis,roi,imagthres,pupilthres]=longaxisdraw_old(mov);
        save(outfile,'longaxis','roi','imagthres','pupilthres');
        disp('saved')
    else
        display(strcat(data_time_stamp_filename,'pupil video not exist'));
    end
%     if exist(strcat(data_time_stamp_filename,'axisroi.mat'),'file')
%         disp('processing pupil after axis')
%         pointdist=@(a,b) sqrt(sum((a-b).^2));
%         warning off;
%         [roiint,areaii,centerx,centery,circularity,framenum]=pupildetect_bw(mov,longaxis,roi,thres,thres2);
%         save(strcat(outfile,'roi.mat'),'roiint');
%         save(strcat(outfile,'pupil.mat'),'areaii','centerx','centery','circularity');
%     else
%         display(strcat(data_time_stamp_filename,'axis roi not exist'));        
%     end
    
end

end