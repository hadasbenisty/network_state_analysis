function main_network_centrality_evaluateion_spont
addpath(genpath('../utils'));
addpath(genpath('D:/utils/affinity/'))
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
animals={'xt' 'xu' 'xs'  'xw' 'xx', 'xz'};%,%
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
similarity_name = { 'pearson_corr'};
signames = {'LSSC' 'Allen'  'Grid4'};% ,};  

for sim_i = 1:length(similarity_name)
    
    for ai = 1:length(animals)
%         eval_weights_and_cent_perday(similarity_name{sim_i}, animals{ai}, statenames);
                eval_weights_and_cent(signames, similarity_name{sim_i}, animals{ai}, statenames);
    end
end
end


function eval_weights_and_cent(signames, simname, animal, statenames)
outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];


mkNewDir(outputfolder);
for s=1:length(signames)
    switch signames{s}
        case 'Allen'
            [~, allen_parcels] = getParcellsByLansAllansAtlas;
            [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
            
            regionLabel.Allen = allen_parcels.regionNum;
            regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
        case 'LSSC'
             [~, allen_parcels] = getParcellsByLansAllansAtlas;
            [parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
            
            regionLabel.Allen = allen_parcels.regionNum;
            regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
            % [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
            [parcels_names.LSSC, ~, ~, regionLabel.LSSC] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
                regionLabel.Allen, animal);
        case 'Grid4'
            load('X:\Hadas\Meso-imaging\lan\xspsych\spt\xs_31_grid4_dfff.mat','par_inds');
            [parcels_names.Grid4, regionLabel.Grid4, finalindex.Grid4, regionLabel.nameslegend, maskByAllen.Grid4, labelsbyallen.Grid4] = getAllenClusteringLabelsGrid(par_inds, 4);
ii = discard_inds;
        
        parcels_names.Grid4=parcels_names.Grid4(ii);
        regionLabel.Grid4=regionLabel.Grid4(ii);
    end
end
load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\',animal,'_spont_data_3states_dfff.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>
disp(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data = eval(statenames{state_i});
    for sig_i = 1:length(signames)
        if 0&exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '.mat']),'file')
        continue;
        end
    
    
        data.(signames{sig_i}) = data.(signames{sig_i})(:, all(~isnan(data.(signames{sig_i}))));
        if any(isnan(data.(signames{sig_i})(:)))
            tt = ~isnan(sum(data.(signames{sig_i})));
            if sum(tt)==0
                disp('nans in dataset')
                
                continue;
            else
                data.(signames{sig_i}) = data.(signames{sig_i})(:,tt);
            end
        end
        
        
        if strcmp(signames{sig_i}, 'Grid4')
            data.(signames{sig_i}) = data.(signames{sig_i})(ii,:);
            thvals = [inf 150:50:400];
        else
            thvals=[inf 5:2:23];
        end
        if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_W.mat']),'file')
            load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_W.mat']),'W_corr')
        else
            W_corr = measure_weights(data.(signames{sig_i}), simname);
            save(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_W.mat']),'W_corr');
        end
        for th=thvals
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}), @process_sim, th);
            
            %         [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_' num2str(th) '.mat']),'W_corr',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
        end
    end
    
end
end




function eval_weights_and_cent_perday(simname, animal, statenames)
outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];


mkNewDir(outputfolder);

[~, allen_parcels] = getParcellsByLansAllansAtlas;
[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;

regionLabel.Allen = allen_parcels.regionNum;
regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
% [roiLabelsbyAllen_gal, regionLabel.Gal, maskByAllen_gal, maskByAllen.Gal] = get_gal_parcels_lables(animal);
[parcels_names.LSSC, ~, ~, regionLabel.LSSC] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
    regionLabel.Allen, animal);

load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\',animal,'_spont_data_3states_dfff.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>
disp(animal);
signames = {'Allen' 'Grid4',};
[~,days2process] = animaltodays(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'Allen.mat']),'file') &&...
            exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'LSSC.mat']),'file')
        continue;
    end
    
    data = eval(statenames{state_i});
    for sig_i = 1:length(signames)
        for dd=1:length(days2process)
            disp([animal ' ' signames{sig_i} ' ' num2str(days2process(dd)) ' ' statenames{state_i}]);
            currdata = data.(signames{sig_i});
            currdata=currdata(:, data.days==days2process(dd));
            if sum(data.days==days2process(dd))==0
                continue;
            end
            if any(currdata(:))
                tt = ~isnan(sum(currdata));
                if sum(tt)==0
                    disp('nans in dataset')
                    
                    continue;
                else
                    currdata = currdata(:,tt);
                end
            end
            if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'file')
                load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'W_corr')
            else
                W_corr = measure_weights(currdata, simname);
            end
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signames{sig_i} '_day' num2str(days2process(dd)) '.mat']),'W_corr',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            
        end
    end
end
end


