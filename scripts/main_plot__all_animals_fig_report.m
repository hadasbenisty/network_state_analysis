function main_plot__all_animals_fig_report

%
addpath('..\functions');
addpath(genpath('../utils'));
addpath(genpath('../meta_data_processing/'));
% procdatapath = 'X:\Hadas\Mesoimaging\crispr\meso_results\ProcessingDirectory_crispr';
[~, ~, finalindex.Allen] = get_allen_meta_parcels;
outputfiggolder = 'X:\Hadas\Mesoimaging\crispr\meso_results\figs_crispr\initial_report';
dffpath = 'X:\Hadas\Meso-imaging\Antara\final_preprocess\alldata';

mkNewDir(outputfiggolder);
animals_db = get_animals_meta_data_by_csv;
for fi = 1:length(animals_db.ex_lut)
    mkNewDir(fullfile(outputfiggolder , ['group' num2str(animals_db.ex_lut(fi))]));
end

spontindic = animals_db.sessionsid_list==find(strcmp(animals_db.sessionsids_lut, 'Spon'));
animals_inds_list=cell(length(animals_db.animal_lut),1);
for i = 1:length(animals_db.animal_lut)
    curranimal = find((animals_db.animal_list ==  (i))&spontindic);
    animals_inds_list{i} = curranimal;
end
LLj=1;
LL=1;figure;
for i = 1:length(animals_inds_list)
    
    for ii = 1:length(animals_inds_list{i})
        k = animals_inds_list{i}(ii);
        if  animals_db.sessionsid_list(k)==find(strcmp(animals_db.sessionsids_lut, 'Air'))
            continue;
        end
        
        
        
        
        datafile_allen = fullfile(dffpath, animals_db.folder_list{k},  'Ca_traces_spt_patch11_Allen_dfff.mat');
        if ~isfile(datafile_allen)
            disp(['No allen '  animals_db.folder_list{k}]);
            continue;
        end
        fulldfof = fullfile(dffpath, animals_db.folder_list{k}, 'Ca_traces_spt_patch11_dfff.mat');
        if ~isfile(fulldfof)
            disp(['No spt '  animals_db.folder_list{k}]);
            continue;
        end
        rawblue = fullfile(dffpath, animals_db.folder_list{k}, 'detrended_blue_dfff_al.mat');
        if ~isfile(rawblue)
            disp(['No raw '  animals_db.folder_list{k}]);
            
            continue;
        end
        
        load(datafile_allen, 'parcels_time_trace');
        
        corr_right = corr(parcels_time_trace(finalindex.Allen,1e2:end-100)');
        corr_left = corr(parcels_time_trace(finalindex.Allen-1,1e2:end-100)');
        
        load('brain_mask','brain_mask');
        load(fulldfof, 'regsig');
        load(rawblue,'sval');
        P=zeros(256);
        
       
        P(brain_mask==1) = sval;
        str = animals_db.animal_lut{animals_db.animal_list(k)};
        str(str == '_') = '-';
        str = [str ' G #' num2str(animals_db.ex_lut(animals_db.ex_list(k)))];
        subplot(5,4,LL);imagesc(P);title(str);colorbar;
        subplot(5,4,LL+1);imagesc(reshape(regsig(:,1e3),[256 256]));title('regressed');colorbar;
        
        subplot(5,4,LL+2);imagesc(corr_right,[0 1]);colorbar;title('Parcels right');
        subplot(5,4,LL+3);imagesc(corr_left,[0 1]);colorbar;title('Parcels left');
        
        set(gcf,'Position',[1          41        1920         963]);
        LL=LL+4;
        if LL>20
            mysave(gcf, fullfile(outputfiggolder, ['summary'  num2str(LLj)]));
            figure;
            LL=1;
            LLj=LLj+1;
        end
        
        
        break;
    end
end
end

