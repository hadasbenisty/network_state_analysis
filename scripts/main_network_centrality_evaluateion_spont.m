function main_network_centrality_evaluateion_spont
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
animals={'xt','xu' 'xs' 'xx','xz','xw'};% 
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
similarity_name = { 'corr', 'fullcorr'};%
for sim_i = 1:length(similarity_name)
    
    for ai = 1:length(animals)
        eval_weights_and_cent(similarity_name{sim_i}, animals{ai}, statenames);
    end
end
end



function eval_weights_and_cent(simname, animal, statenames)
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
signames = {'Allen','LSSC'};
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'Allen.mat']),'file') &&...
            exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,'LSSC.mat']),'file')
        continue;
    end
    
    data = eval(statenames{state_i});
    for sig_i = 1:length(signames)
        data.(signames{sig_i}) = data.(signames{sig_i})(:, all(~isnan(data.Allen)));
        if any(isnan(data.(signames{sig_i})(:)))
            tt = ~isnan(sum(data.(signames{sig_i})));
            if sum(tt)==0
                disp('nans in dataset')
                
                continue;
            else
                data.(signames{sig_i}) = data.(signames{sig_i})(:,tt);
            end
        end
        if exist(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signames{sig_i} '.mat']),'file')
            load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signames{sig_i} '.mat']),'W_corr')
        else
            W_corr = measure_weights_partial(data.(signames{sig_i}), simname);
        end
        [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
        save(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signames{sig_i} '.mat']),'W_corr',...
            'cent_corr_weighted',...
            'cent_corr_notweighted', 'G_corr', 'names_corr');
        
    end
    
end
end





