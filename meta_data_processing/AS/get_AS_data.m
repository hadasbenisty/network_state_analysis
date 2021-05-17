function [imagingData, imaging_time, wheel_speed, pupil_Norm, face_Norm, ...
    wheel_time, pupil_time, vis, airpuff, regionLabel_gal, roiLabelsbyAllen_gal, maskByGal] = get_AS_data(currfolder)
fsspike2=5e3;
[~, ~, finalindex] = get_allen_meta_parcels;
regionLabel_gal = [];
roiLabelsbyAllen_gal = [];
maskByGal=[];


%load imaging time series
load(fullfile(currfolder,'final_dFoF_parcels.mat'),'dFoF_parcells');
imagingData.Allen_green = dFoF_parcells.green;
imagingData.Allen_green = imagingData.Allen_green(finalindex,:);
imagingData.Allen_blue = dFoF_parcells.blue;
imagingData.Allen_blue = imagingData.Allen_blue(finalindex,:);
%% load LSSC
ii=strfind(currfolder, 'DualMice\grab');
jj=strfind(currfolder,'\imaging');
kk=strfind(currfolder,'excitation\');
LSSCfile = dir(fullfile( 'X:\Hadas\Meso-imaging\GRABS_Data\LSSC\', currfolder{1}(ii+8:jj), 'imagingwith575excitation',...
    currfolder{1}(kk+11:end), 'LSCC_dfof.mat'));
if ~isempty(LSSCfile)
    
    
LSSCimaging = load(fullfile(LSSCfile.folder, LSSCfile.name), 'imaging_traces_green');
load(fullfile(LSSCfile.folder, LSSCfile.name),'regionLabel_gal',...
    'roiLabelsbyAllen_gal','maskByGal');
imagingData.LSCC_green = LSSCimaging.imaging_traces_green;
end


% load spike2 and pupil/face data
load(fullfile(currfolder,'smrx_signals.mat'),'channels_data','timestamps','timing')
files=dir(currfolder);
for t=1:length(files)
    if contains(files(t).name,'proc','IgnoreCase',true)
        files=files(t);break;
    end
end
load(fullfile(currfolder,files.name),'proc')
pupil_Norm=proc.output.pupilNorm;
face_Norm=proc.output.facePC1CorrNorm;
wheel_speed = channels_data.wheelspeed;

% get pupil wheel and imaging times
wheel_time = (1:length(wheel_speed))/fsspike2;
imaging_time = timestamps.timaging;
pupil_time = timing.pupilcamstart(1:length(pupil_Norm));

vis = zeros(size(pupil_Norm));
airpuff = zeros(size(pupil_Norm));
if isfield(timing, 'stimstart')
for k=1:length(timing.stimstart)
i = findClosestDouble(pupil_time, timing.stimstart(k));
j = findClosestDouble(pupil_time, timing.stimend(k));
vis(i:j) = 1;
end
end

if isfield(timing, 'airpuffstart')
for k=1:length(timing.airpuffstart)
i = findClosestDouble(pupil_time, timing.airpuffstart(k));
j = findClosestDouble(pupil_time, timing.airpuffstart(k));
airpuff(i:j) = 1;
end
end