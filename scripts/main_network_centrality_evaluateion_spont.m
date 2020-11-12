function main_network_centrality_evaluateion_spont
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../graphs_analysis'));
saveplots = false;
animals={'xs','xx','xz','xw','xt','xu'};
statenames = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
statenames_gal = {'low_pup_q_gal', 'high_pup_q_gal', 'high_pup_l_gal'};

outputfiggolder = 'X:\Lav\ProcessingDirectory\parcor_undirected\';
for ai = 1:length(animals)
    %     eval_weights_and_cent_gal(animals{ai}, saveplots, statenames_gal);
    %     eval_weights_and_cent(animals{ai}, saveplots, statenames);
end
% eval_diff_map(animals, statenames);
% plot_centrality_res_gal(animals, outputfiggolder, statenames_gal);
plot_centrality_res(animals, outputfiggolder, statenames);
end

function eval_diff_map(animals, statenames)
addpath(genpath('../centrality_measures/'));
for ai=1:length(animals)
    animal=animals{ai};
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
    for state_i = 1:length(statenames)
        load(strcat(outputfolder,'network_analysis_corr',statenames{state_i}),'W_corr');
        W_corr = threshold_cor_matrix(W_corr);
        W_corr = W_corr + eye(size(W_corr));
        [M(:, state_i, ai),Q(state_i, ai)]=community_louvain(abs(W_corr));
        P(:, state_i, ai)=participation_coef(abs(W_corr),M(:, state_i, ai),0);
        
        configParams.maxInd=20;
        [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(exp(W_corr), configParams);
        eigenvals(:, state_i, ai) = Lambda(2:end);
        firsteigenvec(:, state_i, ai) = Psi(:,2);
    end
end
M = mean(eigenvals, 3);
S = std(eigenvals, [],3)/sqrt(size(eigenvals,3)-1);
barwitherr(S,M)
end
function plot_centrality_res_gal(animals, outputfiggolder, statenames)
[~, ~, finalindex] = get_allen_meta_parcels;
cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = zeros(256,256,length(animals));
    end
end
for i=1:length(animals)
    animal=animals{i};
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020',animal);
    for state_i = 1:length(statenames)
        load(fullfile(outputfolder,['network_analysis_corr',statenames{state_i}]));
        cent_features = fieldnames(cent_corr_notweighted);
        for cent_i = 1:length(cent_features)
            P = scores_to_heatmap_gal(cent_corr_notweighted.(cent_features{cent_i}), animal);
            spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i})(:, :, i) = P;
            
        end
        
    end
end
mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
for ni = 1:length(cent_features)-1
    %difference maps, not weighted, for each centrality measure
    braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
    parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
    
    diffmask = nanmean(spon_states_notweighted.high_pup_l_gal.(cent_features{ni}),3)-...
        nanmean(spon_states_notweighted.low_pup_q_gal.(cent_features{ni}),3);
    plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup loc minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
    mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_run_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
    diffmask = nanmean(spon_states_notweighted.high_pup_q_gal.(cent_features{ni}),3)-...
        nanmean(spon_states_notweighted.low_pup_q_gal.(cent_features{ni}),3);
    plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup q minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
    mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_pupilhigh_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
end

end

function plot_centrality_res(animals, outputfiggolder, statenames)
[parcels_names] = get_allen_meta_parcels;
cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = [];
        %         spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = [];
    end
end
for i=1:length(animals)
    animal=animals{i};
    outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020',animal);
    for state_i = 1:length(statenames)
        load(fullfile(outputfolder,['network_analysis_corr',statenames{state_i}]));
        cent_features = fieldnames(cent_corr_notweighted);
        for cent_i = 1:length(cent_features)
            spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
                spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}), ...
                cent_corr_notweighted.(cent_features{cent_i}));
            
            %             spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
            %                 spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}), ...
            %                 cent_corr_weighted.(cent_features{cent_i}));
        end
        
    end
end
% mkNewDir(fullfile(outputfiggolder, 'weighted'))
mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
legstr = {'Low Q', 'High Q', 'Loc'};
for ni = 1:length(cent_features)-1
    graph_overlay_allen_2conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        'spon_threestates',cent_features{ni},['2 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr([1 3]));
    
    
    
    %     graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'weighted'), spon_states_weighted.low_pup_q.(cent_features{ni}),...
    %         spon_states_weighted.high_pup_q.(cent_features{ni}), spon_states_weighted.high_pup_l.(cent_features{ni}),...
    %         'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
    %
    graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        spon_states_notweighted.high_pup_q.(cent_features{ni}), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
    
    %difference maps, not weighted, for each centrality measure
    braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
    parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
    
    %high vs low pup
    graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_q.(cent_features{ni}),...
        spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
    
    graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_q.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)']);
    
    %locomotion vs low pup
    graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
    
    graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
        'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)']);
    
    %locomotion vs high pup
    graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
        spon_states_notweighted.high_pup_q.(cent_features{ni}),'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
    
    graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
        spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.high_pup_q.(cent_features{ni}),...
        'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)']);
    
end

end
function eval_weights_and_cent_gal(animal, saveplots, statenames)
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_3states.mat'],...
    'low_pup_q_gal','low_pup_q_t','high_pup_q_gal','high_pup_q_t','high_pup_l_gal','high_pup_l_t','days_to_process'); %#ok<NASGU>


[parcels_names, ~, finalindex] = get_allen_meta_parcels;
[roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal);

[parcels_names_gal, finalindex_gal, maskByGal] = get_gal_meta_parcels_by_allen(parcels_names, finalindex,...
    roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal);



disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'.mat'),'file')
        continue;
    end
    data = eval([statenames{state_i} ]);
    
    data = data(finalindex_gal, :);
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

function eval_weights_and_cent(animal, saveplots, statenames)
outputfolder=fullfile('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\');
mkNewDir(outputfolder);


load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'spon_3states.mat'],...
    'low_pup_q','low_pup_q_t','high_pup_q','high_pup_q_t','high_pup_l','high_pup_l_t','days_to_process'); %#ok<NASGU>
[parcels_names, parcels_region_labels, finalindex] = get_allen_meta_parcels;
disp(animal)
for state_i = 1:length(statenames)
    disp(statenames{state_i})
    if 0&&exist(strcat(outputfolder,'network_analysis_corr',statenames{state_i} ,'.mat'),'file')
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





