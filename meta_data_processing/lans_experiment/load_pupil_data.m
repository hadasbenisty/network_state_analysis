function [pupil_time, pupil_Norm] = load_pupil_data(animal_name, selectedday, pupil_frame, fsspike2, pth)
if ~exist('pth','var')
    pth = 'X:\Older\Lan\Meso-imaging\';
end
pupil_data = fullfile(pth, animal_name, [animal_name '_D' ...
    num2str(selectedday)  '_pupil_clean.mat']);
if ~exist(pupil_data, 'file')
    error('no pupil_data');
end
pupilarea = load(pupil_data,  'areaii');
pupil_Norm=pupilarea.areaii; %get pupil camera timepoints
[timing.pupilstart,timing.pupilend]=squaredetect(pupil_frame,.5);
if size(pupil_Norm,1)==length(timing.pupilstart(1:end))
    timing.pupilstart=timing.pupilstart(1:end);
elseif size(pupil_Norm,1)<length(timing.pupilstart(1:end))
    timing.pupilstart=timing.pupilstart(1:size(pupil_Norm,1));
    timing.pupilend=timing.pupilend(1:size(pupil_Norm,1));
elseif size(pupil_Norm,1)>length(timing.pupilstart(1:end))
    pupil_Norm=pupil_Norm(1:length(timing.pupilstart),:);
end
pupil_time=timing.pupilstart/fsspike2;