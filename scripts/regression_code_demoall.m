
clear;
cd('C:\Users\SuperComp-HigleyLab\Documents\Lav\network_state_analysis\scripts')
addpath(genpath('../pre_processing'));
addpath(genpath('../meta_data_processing'));
addpath(genpath('../network_state_analysis/utils'));
addpath(genpath('../network_state_analysis/functions'));
addpath(genpath('../network_state_analysis/scripts'));

datapath = 'X:\Hadas\Meso-imaging\Antara\data\Antara\AnalyzedData';
imagingpath='X:\Hadas\Meso-imaging\Antara\final_preprocess\alldata';
spike2path = 'X:\CardinLab\Antara\AnalyzedData\';
animals_db = get_animals_meta_data_by_csv;
procdatapath = 'X:\Hadas\Meso-imaging\crispr\meso_results\ProcessingDirectory';
%% regression and ridge regression
%k=4,9,10
typesvals = unique(animals_db.type_list);
%NaN(4,56,1);
for ti = 1:length(typesvals)
    all_results_type=[];
    curtype = animals_db.type_lut{ti};
    animalsinds = find(animals_db.type_list==ti);
    for kj=1:length(animalsinds)
        k=animalsinds(kj);
        if animals_db.isgoodpupil_list(k)==find(strcmp(animals_db.isgoodpupil_lut, 'GOOD'))&&animals_db.isimagingood_list(k)==find(strcmp(animals_db.isimagingood_lut, 'GOOD'))
            %extract_sustained_state(datapath, procdatapath, spike2path, animals_db.folder_list{k});   
            tic
            isair=animals_db.sessionsid_list(k)==1;
            [outputMat_gcamp_regular, outputMat_grabs_ridge] = ...
                returnCVRSquared(animals_db.folder_list{k},k,fullfile(datapath,animals_db.folder_list{k}),fullfile(imagingpath,animals_db.folder_list{k}),fullfile(spike2path,animals_db.folder_list{k}),isair,25);
            toc
            save(fullfile(datapath,animals_db.folder_list{k},'regression_results'),'outputMat_gcamp_regular');
            if ~isempty(outputMat_gcamp_regular)
            all_results_type=cat(3,all_results_type,outputMat_gcamp_regular);
            else
                disp(strcat('EMPTY for',animals_db.folder_list{k}));
            end
        else
            disp(strcat('data isnt good for',animals_db.folder_list{k}));
        end
        clearvars outputMat_grabs_regular isair
    end
    save(fullfile(datapath,strcat(curtype,'all_regression_results')),'all_results_type');
    clearvars all_results_type
end

load('X:\Hadas\Meso-imaging\lan\results\code\Functions\brain_mask.mat')
    parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
outputfiggolder = 'X:\Hadas\Meso-imaging\crispr\meso_results\figs\regression';
for ti = 1:length(typesvals)
    curtype = animals_db.type_lut{ti};
    allpertype=load(fullfile(datapath,strcat(curtype,'all_regression_results')));
    %order: allVars,pupil,face,wheel;
    %all
    km=size(allpertype.all_results_type,3);
    figure;
    regression_heatmap(brain_mask,parcelsallen.parcells_new.indicators,squeeze(allpertype.all_results_type(1,:,:)),strcat('All Vars',curtype,num2str(km),'sessions'),0.25)
   	mysave(gcf, fullfile(outputfiggolder, strcat(curtype,'allvars')));
    %Pupil
    regression_heatmap(brain_mask,parcelsallen.parcells_new.indicators,squeeze(allpertype.all_results_type(2,:,:)),strcat('Pupil',curtype,num2str(km),'sessions'),0.1)
   	mysave(gcf, fullfile(outputfiggolder, strcat(curtype,'Pupil')));    
    regression_heatmap(brain_mask,parcelsallen.parcells_new.indicators,squeeze(allpertype.all_results_type(3,:,:)),strcat('Face',curtype,num2str(km),'sessions'),0.25)
   	mysave(gcf, fullfile(outputfiggolder, strcat(curtype,'Face')));        
    regression_heatmap(brain_mask,parcelsallen.parcells_new.indicators,squeeze(allpertype.all_results_type(4,:,:)),strcat('Wheel',curtype,num2str(km),'sessions'),0.1)
    mysave(gcf, fullfile(outputfiggolder, strcat(curtype,'Wheel')));    
   
    clearvars allpertype
end