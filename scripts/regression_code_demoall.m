
clear;
cd('C:\Users\SuperComp-HigleyLab\Documents\Lav\network_state_analysis\scripts')
addpath(genpath('../pre_processing'));
addpath(genpath('../meta_data_processing'));
addpath(genpath('../../utils'));
addpath('X:\Hadas\Meso-imaging\Antara\results\code\regression')
datapath = 'X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData';
imagingpath='X:\Hadas\Meso-imaging\Antara\final_preprocess\alldata';
spike2path = 'X:\CardinLab\Antara\AnalyzedData\';
animals_db = get_animals_meta_data_by_csv;
procdatapath = 'X:\Hadas\Meso-imaging\crispr\meso_results\ProcessingDirectory';
%% regression and ridge regression
%outputMat_grabs_regular=
%outputMat_grabs_ridge=? size
typesvals = unique(animals_db.type_list);
for ti = 1:length(typesvals)
    curtype = animals_db.type_lut{ti};
    animalsinds = find(animals_db.type_list==ti);
    all_results_type=[];
    for kj=1:length(animalsinds)
        k=animalsinds(kj);
        if animals_db.isgoodpupil_list(k)==find(strcmp(animals_db.isgoodpupil_lut, 'GOOD'))&&animals_db.isimagingood_list(k)==find(strcmp(animals_db.isimagingood_lut, 'GOOD'))
            %extract_sustained_state(datapath, procdatapath, spike2path, animals_db.folder_list{k});
            tic
            isair=animals_db.sessionsid_list(k)==1;
            [outputMat_gcamp_regular, outputMat_grabs_ridge] = ...
                returnCVRSquared(fullfile(datapath,animals_db.folder_list{k}),fullfile(imagingpath,animals_db.folder_list{k}),fullfile(spike2path,animals_db.folder_list{k}),isair,25);
            toc
            save(fullfile(datapath,animals_db.folder_list{k},'regression_results'),'outputMat_gcamp_regular');
            all_results_type=cat(3,all_results_type,outputMat_gcamp_regular);
        end
        clearvars outputMat_grabs_regular isair
    end
    save(fullfile(datapath,strcat(curtype,'all_regression_results')),'all_results_type');
    clearvars all_results_type
end



