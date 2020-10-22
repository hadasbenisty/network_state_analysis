function main_network_centrality_evaluateion_spont
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
saveplots = false;
animals={'xs','xx','xz','xw','xt','xu'};
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
for ai = 1:length(animals)
    eval_weights_and_cent(animals{ai}, saveplots, statenames);
end
plot_centrality_res(animals, outputfiggolder, statenames);
end


function plot_centrality_res(animals, outputfiggolder, statenames)
[parcels_names] = get_allen_meta_parcels;
cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector'};
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = [];
        spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = [];
    end
end
for i=1:length(animals)
    animal=animals{i};
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020',animal);
    for state_i = 1:length(statenames)
        load(fullfile(outputfolder,['network_analysis_corr',statenames{state_i}]));
        cent_features = fieldnames(cent_corr_weighted);
        for cent_i = 1:length(cent_features)
            spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
                spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}), ...
                cent_corr_notweighted.(cent_features{cent_i}));
            
            spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
                spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}), ...
                cent_corr_weighted.(cent_features{cent_i}));
        end
        
    end
end
mkNewDir(fullfile(outputfiggolder, 'weighted'))
mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
legstr = {'Low Q', 'High Q', 'Loc'};
for ni = 1:length(cent_features)
    
    graph_overlay_allen_3conditions(fullfile(outputfiggolder, 'weighted'), spon_states_weighted.low_pup_q.(cent_features{ni}),...
        spon_states_weighted.high_pup_q.(cent_features{ni}), spon_states_weighted.high_pup_l.(cent_features{ni}),...
        'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
    
    graph_overlay_allen_3conditions(fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        spon_states_notweighted.high_pup_q.(cent_features{ni}), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
    
end

end

function eval_weights_and_cent(animal, saveplots, statenames)
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_3states.mat'],...
    'low_pup_q','low_pup_q_t','high_pup_q','high_pup_q_t','high_pup_l','high_pup_l_t','days_to_process'); %#ok<NASGU>
[parcels_names, parcels_region_labels, finalindex] = get_allen_meta_parcels;
disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if exist(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'.mat'),'file')
        continue;
    end
    data = eval(statenames{state_i});
    data = data(finalindex, :);
    data = data(:, all(~isnan(data)));
    if any(isnan(data(:)))
        disp('nans in dataset')
        continue;
    end
    %% measure weights
    W_corr = measure_weights_partial(data, 'corr');
    disp('corweights done')
    % Graph Analysis
    [indic_corr_weighted, indic_corr_notweighted, cent_corr_weighted, cent_corr_notweighted, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names);
    save(strcat(outputfolder,'network_analysis_corr',statenames{state_i}),'W_corr',...
        'indic_corr_weighted','indic_corr_notweighted','cent_corr_weighted',...
        'cent_corr_notweighted', 'G_corr', 'names_corr');
    
    disp('graph analysis saved')
    if saveplots
        %% Visualization
        % plot graphs
        plot_graph(G_corr, cent_corr_weighted, names_corr, parcels_region_labels);
        suptitle('Correlation');
        set(gcf, 'Position',  [150,150, 1000,500]);
        mysave(gcf, strcat(outputfolder,'correlation_graph',statenames{state_i}), 'fig');
        
        figure;
        for k = 1:length(names_corr)
            subplot(2,4,k);
            bar(cent_corr_weighted.(names_corr{k}));
            set(gca, 'XTickLabel', parcels_names);
            set(gca,'XTickLabelRotation',45);
            set(gcf, 'Position',  [150,150, 1000,500]);
            title(names_corr{k});
        end
        suptitle('Correlation');
        mysave(gcf, strcat(outputfolder,'correlation_centrality',statenames{state_i}), 'fig');
        disp('cntrality plotted and saved')
        
        % plot centrality by node population
        figure;
        subplot(2,1,1);
        bar(sum(indic_corr_weighted,2)/size(indic_corr_weighted, 2));set(gca,'XTickLabel',names_corr);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
        title('Correlation');
        
        mysave(gcf, strcat(outputfolder,'centrality_by_node_pop',statenames{state_i}), 'fig');
    end
    
    
end
end





