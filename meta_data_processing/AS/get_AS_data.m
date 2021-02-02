function [imagingData, imaging_time, wheel_speed, pupil_Norm, face_Norm, ...
    wheel_time, pupil_time, regionLabel_gal, roiLabelsbyAllen_gal, maskByGal] = get_AS_data(animalDir, folder, Condition)
fsspike2=5e3;
[~, ~, finalindex] = get_allen_meta_parcels;

folders=dir(animalDir);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
DirFolders= folders(dirFlags);
noDrugFlag=~contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
preDrugFlag=contains({DirFolders.name},'preDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
postDrugFlag=contains({DirFolders.name},'postDrug','IgnoreCase',true)&~contains({DirFolders.name},'Extra','IgnoreCase',true)& ~contains({DirFolders.name},'vis','IgnoreCase',true) ;%ignore visual stim sessions
switch Condition
    case 'NoDrug'
        DirFolders=DirFolders(noDrugFlag);
    case 'PreDrug'
        DirFolders=DirFolders(preDrugFlag);
    case 'PostDrug'
        DirFolders=DirFolders(postDrugFlag);
end

%for each subfolder
currfolder=fullfile(animalDir, DirFolders(folder).name);
%load imaging time series
load(fullfile(currfolder,'final_dFoF_parcels.mat'),'dFoF_parcells');
imagingData.Allen = dFoF_parcells.green;
imagingData.Allen = imagingData.Allen(finalindex,:);
%% load LSSC
ii=strfind(currfolder, 'DualMice\grab');
jj=strfind(currfolder,'\imaging');
LSSCfile = dir(fullfile( 'X:\Hadas\Meso-imaging\GRABS_Data\LSSC\', currfolder(ii+8:jj), 'imagingwith575excitation',...
    DirFolders(folder).name, 'LSCC_dfof.mat'));

LSSCimaging = load(fullfile(LSSCfile.folder, LSSCfile.name), 'imaging_traces_green');
load(fullfile(LSSCfile.folder, LSSCfile.name),'regionLabel_gal',...
    'roiLabelsbyAllen_gal','maskByGal');
imagingData.LSCC = LSSCimaging.imaging_traces_green;



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
