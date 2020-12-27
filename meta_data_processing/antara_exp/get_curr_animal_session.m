function saveoutputpath = get_curr_animal_session(pth2data, animalid, session_i)

%pth2data='/gpfs/gibbs/pi/cardin-higley/Hadas/data/Antara';
T = readtable('../meta_data_processing/antara_exp/Processing_Pipeline_Full.csv');

try
folders = dir(pth2data);
folders=folders(3:end);
animalidname = folders(animalid).name;
pathname = fullfile(pth2data, animalidname);
folders = dir(pathname);
folders=folders(3:end);
saveoutputpath = fullfile(pathname, folders(session_i).name);
catch
   saveoutputpath=[];
end
    