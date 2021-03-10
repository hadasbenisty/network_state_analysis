
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
typesvals = unique(animals_db.type_list);
exvals = unique(animals_db.ex_list);
all_results_type=NaN(length(typesvals),length(exvals));
for ti = 1:length(typesvals) %for each category 
    %curtype = animals_db.type_lut{ti}; %current type 
    for j=1:length(exvals) % for each expression value 
       %cj=animals_db.ex_lut(j); 
       new_inds = find(animals_db.type_list==ti&animals_db.ex_list==j&animals_db.isimagingood_list==1);
       %see how many animals ==that "type" and have good imaging data and
       %meet that "expression" number on the list 
       unique_animals=length(unique(animals_db.animal_list(new_inds)));
       %find unique animals
       %disp(length(find(animals_db.type_list==ti)))
       all_results_type(ti,j)=unique_animals;
       clearvars new_inds cj
    end
end


sexvals = unique(animals_db.sex_list);
all_results_type_sex=NaN(length(typesvals),length(sexvals),length(exvals));
for ti = 1:length(typesvals) %for each category
    %curtype = animals_db.type_lut{ti}; %current type
    for sx=1:length(sexvals) %for one of each sex 
        %cursex=animals_db.sex_lut(sx);
        for j=1:length(exvals) %for one of each expression type
            %cj=animals_db.ex_lut(j);
            new_inds = find(animals_db.type_list==ti&animals_db.sex_list==sx&animals_db.ex_list==j&animals_db.isimagingood_list);
            % get new indices
            unique_animals=length(unique(animals_db.animal_list(new_inds))); %get indices
            all_results_type_sex(ti,sx,j)=unique_animals; %for each category, store # of unique animals
            clearvars new_inds cj
        end
    end
end
female=squeeze(all_results_type_sex(:,2,:))
male=squeeze(all_results_type_sex(:,3,:))