function main_network_centrality_evaluateion_spont_2p_lan
addpath(genpath('../utils'));
addpath(genpath('D:/utils/affinity/'))
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animals = T.AnimalID;
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
similarity_name = {'pearson_corr',  };%'corr',,  'L2' 'fullcorr' 'cov''partial_corr'
signames = {'cells'  };% ,'LSSC'};'Allen'

for sim_i = 1:length(similarity_name)
    
    for ai = 1:length(animals)
        eval_weights_and_cent_perday(signames, similarity_name{sim_i}, animals{ai}, statenames);
        %         eval_weights_and_cent(signames, similarity_name{sim_i}, animals{ai}, statenames);
    end
end
end

function eval_weights_and_cent_perday(signames, simname, animal, statenames)
outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];


mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\',animal,'_spont_data_3states_dfff.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>
disp(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data = eval(statenames{state_i});
    dayslist = unique(data.days);
    
    for sig_i = 1:length(signames)
        for day_i = 1:length(dayslist)
            if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_D' num2str(dayslist(day_i)) '_' num2str(700) '.mat']),'file')
                continue;
            end
            if isempty(data.(signames{sig_i}) )
                continue;
            end
            daysinds = data.days == dayslist(day_i);
            currdata = data.(signames{sig_i})(:,daysinds);
            currdata = currdata(:, all(~isnan(currdata)));
            if any(isnan(currdata(:)))
                tt = ~isnan(sum(currdata));
                if sum(tt)==0
                    disp('nans in dataset')
                    
                    continue;
                else
                    currdata = currdata(:,tt);
                end
            end
            
            
            if 0& exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_D' num2str(dayslist(day_i)) '.mat']),'file')
                load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_D' num2str(dayslist(day_i)) '.mat']),'W_corr')
            else
                W_corr = measure_weights(currdata, simname);
            end
            ththird = round(size(W_corr,1));
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] =...
                graph_analysis_afterclust(W_corr, 1:size(W_corr), [], @process_sim, ththird);
            
            %         [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_D' num2str(dayslist(day_i)) '_ththird.mat']),'W_corr',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
            thT = size(W_corr,1);
%             for th=3:2:thT
%                 [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] =...
%                     graph_analysis_afterclust(W_corr, 1:size(W_corr), [], @process_sim, th);
%                 
%                 %         [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
%                 save(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_D' num2str(dayslist(day_i)) '_' num2str(th) '.mat']),'W_corr',...
%                     'cent_corr_weighted',...
%                     'cent_corr_notweighted', 'G_corr', 'names_corr');
%             end
        end
    end
end
end


function eval_weights_and_cent(signames, simname, animal, statenames)
outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\network_centrality_' simname];


mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\',animal,'_spont_data_3states_dfff.mat'],...
    'low_pup_q','high_pup_q','high_pup_l'); %#ok<NASGU>
disp(animal);
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    
    data = eval(statenames{state_i});
    for sig_i = 1:length(signames)
        if exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_' num2str(3) '.mat']),'file')
            continue;
        end
        if isempty(data.(signames{sig_i}) )
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
        
        
        
        if 0& exist(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '.mat']),'file')
            load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '.mat']),'W_corr')
        else
            W_corr = measure_weights(data.(signames{sig_i}), simname);
        end
        thT = size(W_corr,1);
        for th=3:2:thT
            [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] =...
                graph_analysis_afterclust(W_corr, 1:size(W_corr), [], @process_sim, th);
            
            %         [cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names.(signames{sig_i}), regionLabel.(signames{sig_i}));
            save(fullfile(outputfolder,[animal '_',statenames{state_i} ,signames{sig_i} '_' num2str(th) '.mat']),'W_corr',...
                'cent_corr_weighted',...
                'cent_corr_notweighted', 'G_corr', 'names_corr');
        end
    end
    
end
end





