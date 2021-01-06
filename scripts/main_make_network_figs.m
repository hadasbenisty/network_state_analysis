function main_make_network_figs
addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
animals={'xt','xu' 'xs', 'xx','xz','xw'};%,
stateslabels = { 'low_pup_q', 'high_pup_q', 'high_pup_l'};
cent_features = {  'eigenvector' 'degree' 'closeness' 'participation' , 'diffmap',  'betweenness' 'pagerank', 'second_eigval'};%};%'eigenvector'
similarity_name = {'partial_corr'   };%'pearson_corr', 'corr',,  'fullcorr' 'cov''partial_corr'
doover=false;
%% Fig 4 - network
% plot_centrality_res_per_day(animals, outputfiggolder, stateslabels);
for sim_i = 1:length(similarity_name)
    outputfiggolder = ['X:\Hadas\Meso-imaging\lan\meso_results\figs\network_centrality_' similarity_name{sim_i}];
    mkNewDir(outputfiggolder)
    %     plot_similarity_res(similarity_name{sim_i}, animals, outputfiggolder, stateslabels);
    
    plot_centrality_res(cent_features, similarity_name{sim_i}, animals, outputfiggolder, stateslabels, doover);
    
end
end


function [trial_states, spon_states] = load_similarity_matrices(signame, outputfolder, animals, statenames, isweigtedstr)

for i=1:length(animals)
    animal=animals{i};
    for state_i = 1:length(statenames)
        % trial correct
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame]), ...
            'W_corr_cor');
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'trials_incorrect_' signame]), ...
            'W_corr_inc');
        load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signame]), ...
            'W_corr');
        switch signame
            case 'Allen'
                
                
                spon_states(:,:,state_i, i) = process_sim(W_corr); %#ok<AGROW>
                trial_states.correct(:,:,state_i, i) = process_sim(W_corr_cor);
                trial_states.incorrect(:,:,state_i, i) = process_sim(W_corr_inc);
            case 'LSSC'
                spon_states.(statenames{state_i}){i}  = process_sim(W_corr);
                trial_states.(statenames{state_i}).correct{i} = process_sim(W_corr_cor);
                trial_states.(statenames{state_i}).incorrect{i} = process_sim(W_corr_inc);
        end
        
    end
end


end

function [trial_states, spon_states, spont_heatmap, trial_heatmap]  = load_centrality_results_per_days_allen(cent_features, signame, outputfolder, animals, statenames, isweigtedstr)
for state_i = 1:length(statenames)
    for cent_i = 1:length(cent_features)
        spon_states.(statenames{state_i}).(cent_features{cent_i}) = nan(23, 40, 6);
        trial_states.(statenames{state_i}).(cent_features{cent_i}).correct = nan(23, 40, 6);
        trial_states.(statenames{state_i}).(cent_features{cent_i}).incorrect =nan(23, 40, 6);
        spont_heatmap.(statenames{state_i}).(cent_features{cent_i}) = zeros(256, 256, length(animals));
        trial_heatmap.(statenames{state_i}).(cent_features{cent_i}).correct = zeros(256, 256, length(animals));
        trial_heatmap.(statenames{state_i}).(cent_features{cent_i}).incorrect = zeros(256, 256, length(animals));
        
    end
end


for i=1:length(animals)
    animal=animals{i};
    [~,days_list] = animaltodays(animal);
    for day_i = 1:length(days_list)
        
        for state_i = 1:length(statenames)
            % trial correct
            filename = fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            trial_states.(statenames{state_i}).(cent_features{cent_i}).correct(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                end
            end
            % trial incorrect
            filename = fullfile(outputfolder,[animal '_',statenames{state_i} ,'trials_incorrect_' signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            trial_states.(statenames{state_i}).(cent_features{cent_i}).incorrect(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                end
                
            end
            % spont
            filename = fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} ,signame '_day' num2str(days_list(day_i)) '.mat']);
            if exist(filename, 'file')
                load(filename, ['cent_corr_' isweigtedstr ]);
                centvals = eval(['cent_corr_' isweigtedstr ]);
                for cent_i = 1:length(cent_features)
                    switch signame
                        case 'Allen'
                            spon_states.(statenames{state_i}).(cent_features{cent_i})(:, days_list(day_i), i) = centvals.(cent_features{cent_i});
                    end
                    
                end
            end
        end
    end
end
end

function [correct_states, incorrect_states, spon_states, correct_heatmap, incorrect_heatmap, spont_heatmap]  = load_centrality_results(cent_features, signame, outputfolder, animals, statenames, isweigtedstr, th)
switch signame
    case 'Allen'
        
        spon_states = nan(23, length(animals), length(statenames), length(cent_features));
        correct_states = nan(23, length(animals), length(statenames), length(cent_features));
        incorrect_states = nan(23, length(animals), length(statenames), length(cent_features));
    case 'LSSC'
        for ai=1:length(animals)
            q=load(['X:\Hadas\Meso-imaging\lan\' animals{ai} 'psych\spt\gal\' animals{ai} '_finalparcellation_Gal.mat'], 'final_par_mask');
            mask=q.final_par_mask;
            mask(:,1:128,:) = 0;
            N=0;
            for n=1:size(mask,3)
                if any(any(mask(:,:,n)>0))
                    N=N+1;
                end
            end
            spon_states{ai} = nan(N,  length(statenames), length(cent_features));
            correct_states{ai} = nan(N,  length(statenames), length(cent_features));
            incorrect_states{ai} = nan(N,  length(statenames), length(cent_features));
        end
    case 'Grid4'
        ii=discard_inds;
        spon_states = nan(length(ii), length(animals), length(statenames), length(cent_features));
        correct_states = nan(length(ii), length(animals), length(statenames), length(cent_features));
        incorrect_states = nan(length(ii), length(animals), length(statenames), length(cent_features));
        
end
spont_heatmap = nan(256,256, length(animals),length(statenames), length(cent_features));
correct_heatmap = nan(256,256, length(animals),length(statenames), length(cent_features));
incorrect_heatmap = nan(256,256, length(animals),length(statenames), length(cent_features));






for i=1:length(animals)
    animal=animals{i};
    for state_i = 1:length(statenames)
        % trial correct
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_correct_' signame '_' num2str(th)]), ...
            ['cent_corr_' isweigtedstr ]);
        centvals = eval(['cent_corr_' isweigtedstr ]);
        xx=[];
        for cent_i = 1:length(cent_features)
            xx(:, cent_i) = centvals.(cent_features{cent_i});
            if strcmp(signame,'Allen')||strcmp(signame,'Grid4')
                correct_states(:, i, state_i, cent_i) = centvals.(cent_features{cent_i});
            else
                correct_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        Pvec = scores_to_heatmap(xx, 0,signame, animal);
        for cent_i = 1:length(cent_features)
            correct_heatmap(:,:,i, state_i,cent_i)=Pvec(:,:,cent_i);
        end
        % trial incorrect
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,'_trials_incorrect_' signame '_' num2str(th)]), ...
            ['cent_corr_' isweigtedstr ]);
        centvals = eval(['cent_corr_' isweigtedstr ]);
        xx=[];
        for cent_i = 1:length(cent_features)
            xx(:, cent_i) = centvals.(cent_features{cent_i});
            if strcmp(signame,'Allen')||strcmp(signame,'Grid4')
                incorrect_states(:, i, state_i, cent_i) = centvals.(cent_features{cent_i});
            else
                incorrect_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        Pvec = scores_to_heatmap(xx, 0,signame, animal);
        for cent_i = 1:length(cent_features)
            incorrect_heatmap(:,:,i, state_i,cent_i)=Pvec(:,:,cent_i);
        end
        
        
        
        % spont
        load(fullfile(outputfolder,[animal '_',statenames{state_i} ,signame '_' num2str(th)]), ...
            ['cent_corr_' isweigtedstr ] );
        centvals = eval(['cent_corr_' isweigtedstr ]);
        
        xx=[];
        for cent_i = 1:length(cent_features)
            xx(:, cent_i) = centvals.(cent_features{cent_i});
            if strcmp(signame,'Allen')||strcmp(signame,'Grid4')
                spon_states(:, i, state_i, cent_i) = centvals.(cent_features{cent_i});
            else
                spon_states{i}(:, state_i, cent_i) = centvals.(cent_features{cent_i});
            end
        end
        Pvec = scores_to_heatmap(xx, 0,signame, animal);
        for cent_i = 1:length(cent_features)
            spont_heatmap(:,:,i, state_i,cent_i)=Pvec(:,:,cent_i);
        end
        
    end
    
end
end
function plot_centrality_res_per_day(animals, outputfiggolder, statenames)

outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';

cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation'};
signals_names = {'Allen', 'LSSC'};
isweigtedstr = {'notweighted' 'weighted'};
%
for sig_i = 1:length(signals_names)
    for isweigted = 1:2
        sumfile = fullfile(outputfolder, ['summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'per_days.mat']);
        if exist(sumfile, 'file')
            load(sumfile, 'trial_states', 'spon_states', 'spont_heatmap', 'trial_heatmap');
        else
            if sig_i==1
                [trial_states, spon_states, spont_heatmap, trial_heatmap] = load_centrality_results_per_days_allen(cent_features, signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted});
            else
                error('')
            end
            save(sumfile, 'trial_states', 'spon_states', 'spont_heatmap', 'trial_heatmap');
        end
        
        legstr = {'Low Q', 'High Q', 'Loc'};
        if sig_i==1
            parcels_names = get_allen_meta_parcels;
            N = size(trial_states.low_pup_q.eigenvector.correct,2);
            for l=1:length(cent_features)
                centstr = cent_features{l};
                for k=1:length(statenames)
                    M(:, 1, k) = nanmean(trial_states.(statenames{k}).(centstr).correct,2);
                    M(:, 2, k) = nanmean(trial_states.(statenames{k}).(centstr).incorrect,2);
                    S(:, 1, k) = nanstd(trial_states.(statenames{k}).(centstr).correct,[],2)/sqrt(N-1);
                    S(:, 2, k) = nanstd(trial_states.(statenames{k}).(centstr).incorrect,[],2)/sqrt(N-1);
                end
                plot_correct_incorrect_per_state_per_parcels(M, S, parcels_names, statenames)
                suptitle(cent_features{l});
                % mysave(gcf, fullfile(outputfiggolder,['trial_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i}]));
                
            end
            for ni = 1:length(cent_features)-1
                graph_overlay_allen_3conditions('', '', spon_states.low_pup_q.(cent_features{ni}),...
                    spon_states.high_pup_q.(cent_features{ni}), spon_states.high_pup_l.(cent_features{ni}),...
                    '',cent_features{ni},['3 states ' cent_features{ni} ' Centrality '], parcels_names,length(animals), legstr);
                % mysave(gcf, fullfile(outputfiggolder,['spont_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i}]));
                
            end
        end
        
        %% heatmaps
        % spont
        for ni = 5% 1:length(cent_features)-1
            figure;
            for state_i = 1:length(statenames)
                subplot(3,3,state_i)
                mask = mean(spont_heatmap.(statenames{state_i}).(cent_features{ni})(:,:,[1 2 3 5 6]),3);
                plot_vals_heatmap(mask, 'Centrality',...
                    [],  0, 0.008, 1,colormap(redblue))
                title([statenames{state_i} ' Spont']);
                mask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).correct(:,:,[1 2 3 5 6]),3);
                subplot(3,3,3+state_i)
                plot_vals_heatmap(mask, 'Centrality',...
                    [],  0, 0.008, 1,colormap(redblue));
                title([statenames{state_i} ' Correct']);
                mask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).incorrect(:,:,[1 2 3 5 6]),3);
                subplot(3,3,6+state_i)
                plot_vals_heatmap(mask, 'Centrality',...
                    [],  0, 0.008, 1,colormap(redblue))
                title([statenames{state_i} ' Incorrect']);
            end
            suptitle(cent_features{ni});
            % mysave(gcf, fullfile(outputfiggolder,[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
            
            
            diffmask = mean(spont_heatmap.high_pup_l.(cent_features{ni})-spont_heatmap.low_pup_q.(cent_features{ni}),3);
            figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
            plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                [],  -L, L, 1,colormap(redblue))
            title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
            % mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
            % trial
            for state_i=1:length(statenames)
                diffmask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).correct-trial_heatmap.(statenames{state_i}).(cent_features{ni}).incorrect,3);
                figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
                plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                    [],  -L, L, 1,colormap(redblue))
                title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
                
                % mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
            end
        end
    end
end
end

% function plot_similarity_res(simname, animals, outputfiggolder, statenames)
%
% outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];
%
% signals_names = {'Allen', 'LSSC'};
% isweigtedstr = {'notweighted' 'weighted'};
%
% for sig_i = 1:length(signals_names)
%     for isweigted = 1:2
%
%
%         [trial_states, spon_states] = load_similarity_matrices(signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted});
%
%
%         legstr = {'Low Q', 'High Q', 'Loc'};
%         if sig_i==1
%             [parcels_names,parcels_region_labels] = get_allen_meta_parcels;
%             figure;
%             for k=1:length(statenames)
%             G = graph(mean(spon_states(:,:,k,:),4),parcels_names);
%             subplot(3,3,k);p=plot(G,'NodeCData', parcels_region_labels);colormap(jet);
%             title(statenames{k});
%
%             G = graph(mean(trial_states.correct(:,:,k,:),4),parcels_names);
%             subplot(3,3,k+3);p=plot(G,'NodeCData', parcels_region_labels);colormap(jet);
%             title('Correct');
%
%             G = graph(mean(trial_states.incorrect(:,:,k,:),4),parcels_names);
%             subplot(3,3,k+6);p=plot(G,'NodeCData', parcels_region_labels);colormap(jet);
%             title('Incorrect');
%             end
%
%             N = size(trial_states.correct,2);
%             figure;
%             for k=1:length(statenames)
%                 subplot(3,3,k);
%                 imagesc(mean(spon_states(:,:,k,:),4));
%                 set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
%                 set(gca,'YTick',1:23);set(gca,'YTickLabel',parcels_names)
%                 title(statenames{k});
%                 subplot(3,3,3+k);
%                 imagesc(mean(trial_states.correct(:,:,k,:),4));
%                 set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
%                 set(gca,'YTick',1:23);set(gca,'YTickLabel',parcels_names)
%                 title('Correct');
%                 subplot(3,3,6+k);
%                 imagesc(mean(trial_states.incorrect(:,:,k,:),4));
%                 set(gca,'XTick',1:23);set(gca,'XTickLabel',parcels_names)
%                 set(gca,'YTick',1:23);set(gca,'YTickLabel',parcels_names);
%                 title('Incorrect');
%             end
%             mysave(gcf, fullfile(outputfiggolder,['trial_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i}]));
%
%         end
%         for ni = 1:length(cent_features)-1
%             graph_overlay_allen_3conditions('', '', spon_states.low_pup_q.(cent_features{ni}),...
%                 spon_states.high_pup_q.(cent_features{ni}), spon_states.high_pup_l.(cent_features{ni}),...
%                 '',cent_features{ni},['3 states ' cent_features{ni} ' Centrality '], parcels_names,length(animals), legstr);
%             mysave(gcf, fullfile(outputfiggolder,['spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i}]));
%
%         end
%     end
%
%     %% heatmaps
%     % spont
%     for ni = 1:length(cent_features)-1
%         figure;
%         mask = mean(spont_heatmap.(statenames{1}).(cent_features{ni}),3);
%         LM =  max(mask(:));
%         Lm =  min(mask(:));
%         for state_i = 1:length(statenames)
%             subplot(3,3,state_i);
%             mask = mean(spont_heatmap.(statenames{state_i}).(cent_features{ni}),3);
%             LM =  max(mask(:));
%             Lm =  min(mask(:));
%             plot_vals_heatmap(mask, 'Centrality',...
%                 [],  Lm, LM, 1,colormap(redblue))
%             title([statenames{state_i} ' Spont']);
%             mask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).correct,3);
%             LM =  max(mask(:));
%             Lm =  min(mask(:));
%             subplot(3,3,3+state_i)
%             plot_vals_heatmap(mask, 'Centrality',...
%                 [],  Lm, LM, 1,colormap(redblue));
%             title([statenames{state_i} ' Correct']);
%             mask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).incorrect,3);
%             LM =  max(mask(:));
%             Lm =  min(mask(:));
%             subplot(3,3,6+state_i)
%             plot_vals_heatmap(mask, 'Centrality',...
%                 [],  Lm, LM, 1,colormap(redblue))
%             title([statenames{state_i} ' Incorrect']);
%         end
%         suptitle(cent_features{ni});
%         mysave(gcf, fullfile(outputfiggolder,[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
%
%
%         diffmask = mean(spont_heatmap.high_pup_l.(cent_features{ni})-spont_heatmap.low_pup_q.(cent_features{ni}),3);
%         figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%         plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%             [],  -L, L, 1,colormap(redblue))
%         title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
%         mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
%         % trial
%         for state_i=1:length(statenames)
%             diffmask = mean(trial_heatmap.(statenames{state_i}).(cent_features{ni}).correct-trial_heatmap.(statenames{state_i}).(cent_features{ni}).incorrect,3);
%             figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%             plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%                 [],  -L, L, 1,colormap(redblue))
%             title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
%
%             mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i}]));
%         end
%     end
%     M=[];S=[];
%     for state_i=1:length(statenames)
%         M(state_i) = nanmean(spon_states.(statenames{state_i}).second_eigval);
%         S(state_i) = nanstd(spon_states.(statenames{state_i}).second_eigval)/sqrt(length(spon_states.(statenames{state_i}).second_eigval)-1);
%     end
%     figure;subplot(3,1,1);plot_3_bars(M,S,statenames);
%     ylim([0 .5]);
%     M=[];S=[];
%     for state_i=1:length(statenames)
%
%         M(state_i,1) = nanmean(trial_states.(statenames{state_i}).second_eigval.correct);
%         S(state_i,1) = nanstd(trial_states.(statenames{state_i}).second_eigval.correct)/sqrt(length(trial_states.(statenames{state_i}).second_eigval.incorrect)-1);
%         M(state_i,2) = nanmean(trial_states.(statenames{state_i}).second_eigval.incorrect);
%         S(state_i,2) = nanstd(trial_states.(statenames{state_i}).second_eigval.incorrect)/sqrt(length(trial_states.(statenames{state_i}).second_eigval.incorrect)-1);
%
%     end
%
%     subplot(3,1,2);
%     plot_3_bars(M(:,1),S(:,1),statenames);ylim([0 0.5]);title('c');
%     subplot(3,1,3);plot_3_bars(M(:,2),S(:,2),statenames);ylim([0 0.5]);title('i');
%     mysave(gcf, fullfile(outputfiggolder,['second_eigval_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i}]));
%
% end
% close all;
%
% end

function plot_centrality_res(cent_features, simname, animals, outputfiggolder, statenames, doover)

outputfolder=['X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality_' simname];
files = dir(['X:\Hadas\Meso-imaging\lan\' animals{1} 'psych\spt\' animals{1} '*_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

signals_names = {  'Allen'};%'LSSC' 'Grid4'

isweigtedstr = { 'weighted'};%'notweighted'
for l = 1:length(cent_features)
    mkNewDir(fullfile(outputfiggolder,cent_features{l}));
end

for sig_i = 1:length(signals_names)
  
        
    switch signals_names{sig_i}
        case 'Allen'
            [parcels_names, parcels_region_labels, final_index, region_lut, grid_map_final_index, labelsbyallen] = get_allen_meta_parcels;
            visualinds = 1;
            somatoinds = 14;
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
             thT=7:2:23;
        case 'Grid4'
            [parcels_names, parcels_region_labels, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
            visualinds = find(parcels_region_labels==1);
            somatoinds = find(parcels_region_labels==6);
            thT=200;%50:50:400;
    end
    for isweigted = 1:length(isweigtedstr)
        for th=thT
            sumfile = fullfile(outputfolder, ['summary_centrality_', signals_names{sig_i}, '_' isweigtedstr{isweigted} 'th' num2str(th) '.mat']);
            if  ~doover&&exist(sumfile, 'file')
                load(sumfile);
            else
                [correct_states, incorrect_states, spon_states, correct_heatmap,...
                    incorrect_heatmap, spont_heatmap] = load_centrality_results(cent_features, signals_names{sig_i}, outputfolder, animals, statenames, isweigtedstr{isweigted}, th);
                save(sumfile, 'correct_states', 'incorrect_states', 'spon_states',...
                    'correct_heatmap', 'incorrect_heatmap', 'spont_heatmap');
            end
            diffmapind = find(strcmp(cent_features, 'diffmap'));
            if ~isempty(diffmapind)
                for j=1:length(statenames)
                    if sig_i==1
                        spon_states(:, :, j, diffmapind) = sign(diffmap2clusters(spon_states(:, :, j, diffmapind)));
                        correct_states(:, :, j, diffmapind) = sign(diffmap2clusters(correct_states(:, :, j, diffmapind)));
                        incorrect_states(:, :, j, diffmapind) = sign(diffmap2clusters(incorrect_states(:, :, j, diffmapind)));
                    end
                    
                    spont_heatmap(:,:,:,j,diffmapind) = sign(diffmap2clusters(spont_heatmap(:,:,:,j,diffmapind)));
                    correct_heatmap(:,:,:,j,diffmapind) = sign(diffmap2clusters(correct_heatmap(:,:,:,j,diffmapind)));
                    incorrect_heatmap(:,:,:,j,diffmapind) = sign(diffmap2clusters(incorrect_heatmap(:,:,:,j,diffmapind)));
                end
            end
            
            legstr = {'Low Q', 'High Q', 'Loc'};
            if ~strcmp(signals_names{sig_i}, 'LSSC')
                for l=find(~strcmp(cent_features, 'second_eigval'))
                    diffspont = spon_states(:, :, 3, l) - spon_states(:,:,1,l);
                    difftrial = correct_states(:,:,:,l) - incorrect_states(:,:,:,l);
                    figure;
                    for state_i = 1:length(statenames)
                        subplot(1,3,state_i);
                        x=mean(difftrial([visualinds somatoinds], :, state_i),2);
                        y=mean(diffspont([visualinds somatoinds], :),2);
                        lmodel = fitlm(x(:), y(:));
                        
                        hold all;
                        for ri = 1:7
                            nds = find(parcels_region_labels==ri);
                            x=mean(difftrial(nds, :, state_i),2);
                        y=mean(diffspont(nds, :),2);
                        plot(x,y,'.', 'MarkerSize',12);
                        end
                        p=plot(lmodel);
                        p(1).Visible='off';
                        legend(region_lut)
                        ylabel([legstr{3} ' minus ' legstr{1} ' spont']);
                        xlabel(['corr-inc on ' legstr{state_i}]);
                       
                        title(sprintf('%s R2=%2.2f', legstr{state_i},  lmodel.Rsquared.Ordinary));
%                         leg = get(gca,'Legend');leg.Visible='off';
                        
                    end
                    suptitle(['Trend of diffs for ' cent_features{l} ' visual and somato'])
                    set(gcf, 'Position',[176         177        1598         627]);
                    mysave(gcf, fullfile(outputfiggolder,cent_features{l}, ['diffs_',cent_features{l},'_'  isweigtedstr{isweigted}  '_' signals_names{sig_i} '_th' num2str(th)]));

                end
            end
            if strcmp(signals_names{sig_i},'Allen')
                parcels_names = get_allen_meta_parcels;
                N = size(correct_states,2);
                Mc = squeeze(nanmean(correct_states,2));
                Mi = squeeze(nanmean(incorrect_states,2));
                Sc = squeeze(nanstd(correct_states,[],2));
                Si = squeeze(nanstd(incorrect_states,[],2));
                for l=find(~strcmp(cent_features, 'second_eigval'))
                    M=[];S=[];
                    for k=1:length(statenames)
                        M(:,1,k) = Mc(:,k,l);
                        M(:, 2, k) = Mi(:,k,l);
                        S(:, 1, k) = Sc(:,k,l)/sqrt(N-1);
                        S(:, 2, k) = Si(:,k,l)/sqrt(N-1);
                    end
                    plot_correct_incorrect_per_state_per_parcels(M, S, parcels_names, statenames)
                    suptitle(cent_features{l});
                    mysave(gcf, fullfile(outputfiggolder,cent_features{l},['trial_',cent_features{l},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
                    
                end
                for ni = find(~strcmp(cent_features, 'second_eigval'))
                    graph_overlay_allen_3conditions('', '', spon_states(:,:,1,ni),...
                        spon_states(:,:,2,ni), spon_states(:,:,3,ni),...
                        '',cent_features{ni},['3 states ' cent_features{ni} ' Centrality '], parcels_names,length(animals), legstr);
                    mysave(gcf, fullfile(outputfiggolder,cent_features{ni}, ['spont_',cent_features{ni},'_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
                    
                end
            end
            
            %% heatmaps
            
            maskscentmeanspont = squeeze(mean(spont_heatmap,3));
            maskscentmeanc = squeeze(mean(correct_heatmap,3));
            maskscentmeani = squeeze(mean(incorrect_heatmap,3));
            for ni = find(~strcmp(cent_features, 'second_eigval'))
                % spont
                %             plot_cent_per_animal(spont_heatmap(:,:,:,:,ni))
                %             suptitle([cent_features{ni} ' Spont']);
                %             mysave(gcf, fullfile(outputfiggolder,cent_features{ni},[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} 'peranimal_spont']));
                % Correct
                %             plot_cent_per_animal(correct_heatmap(:,:,:,:,ni))
                %             suptitle([cent_features{ni} ' Correct']);
                %             mysave(gcf, fullfile(outputfiggolder,cent_features{ni},[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} 'peranimal_correct']));
                %             %Incorrect
                %             plot_cent_per_animal(incorrect_heatmap(:,:,:,:,ni))
                %             suptitle([cent_features{ni} ' Incorrect']);
                %             mysave(gcf, fullfile(outputfiggolder,cent_features{ni},[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} 'peranimal_incorrect']));
                
                
                
                figure;
                v=[reshape(maskscentmeanspont(:,:,:,ni),1,[]) ...
                    reshape(maskscentmeanc(:,:,:,ni),1,[]) reshape(maskscentmeani(:,:,:,ni),1,[]) ];
                L = quantile(v,[.1 .9]);
                if L(1)==L(2)
                    continue;
                end
                for state_i = 1:length(statenames)
                    % spont
                    subplot(4,4,state_i);
                    
                    plot_vals_heatmap(maskscentmeanspont(:,:,state_i,ni), 'Centrality',...
                        [],  L(1), L(2), 1,colormap(redblue))
                    title([statenames{state_i} ' Spont']);
                    
                    subplot(4,4,4+state_i)
                    plot_vals_heatmap(maskscentmeanc(:,:,state_i,ni), 'Centrality',...
                        [],  L(1), L(2), 1,colormap(redblue));
                    title([statenames{state_i} ' Correct']);
                    
                    subplot(4,4,8+state_i)
                    plot_vals_heatmap(maskscentmeani(:,:,state_i,ni), 'Centrality',...
                        [],  L(1), L(2), 1,colormap(redblue))
                    title([statenames{state_i} ' Incorrect']);
                end
                diffmask = maskscentmeanspont(:,:,3,ni)-maskscentmeanspont(:,:,1,ni);
                L = quantile(diffmask(:),[.1 .9]);
                subplot(4,4,4);
                plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                    [],  -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                title('3 - 1 spont');
                
                diffmask = maskscentmeanc(:,:,3,ni)-maskscentmeanc(:,:,1,ni);
                subplot(4,4,8);
                plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                    [],  -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                title('3 - 1 correct');
                
                diffmask = maskscentmeani(:,:,3,ni)-maskscentmeani(:,:,1,ni);
                subplot(4,4,12);
                plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
                    [],  -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                title('3 - 1 incorrect');
                
                diffmasks = maskscentmeanc(:,:,:,ni)-maskscentmeani(:,:,:,ni);
                % L = quantile(diffmasks(:),[.1 .9]);
                for state_i=1:length(statenames)
                    subplot(4,4,state_i+12)
                    plot_vals_heatmap(diffmasks(:,:,state_i), 'Difference in Node Centrality',...
                        [],   -abs(L(1)), abs(L(1)), 1,colormap(redblue))
                    title(['Cor-inc on ' statenames{state_i}  ]);
                    
                    %mysave(gcf, fullfile(outputfiggolder,cent_features{ni},['trial_state' num2str(state_i) '_',cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} '_th' num2str(th)]));
                end
                
                suptitle(cent_features{ni});
                mysave(gcf, fullfile(outputfiggolder,cent_features{ni},[cent_features{ni},'_' isweigtedstr{isweigted} '_heatmap_' signals_names{sig_i} '_th' num2str(th)]));
                
                
                
                
                
                
                
            end
            secondegvali = find(strcmp(cent_features, 'second_eigval'));
            if ~isempty(secondegvali)
                if ~strcmp(signals_names{sig_i},'LSSC')
                    M = squeeze(mean(spon_states(1,:,:,secondegvali),2));
                    S = squeeze(std(spon_states(1,:,:,secondegvali),[],2))/sqrt(size(spon_states,2)-1);
                else
                    v=[];
                    for ai=1:size(spon_states,2)
                        v(:,ai) = squeeze(spon_states{ai}(1,:,secondegvali));
                    end
                    M=mean(v,2);
                    S=std(v,[],2);
                end
                figure;subplot(3,1,1);plot_3_bars(M,S,statenames);
                ylim([0 .2]);
                if ~strcmp(signals_names{sig_i},'LSSC')
                    M(:,1) = squeeze(mean(correct_states(1,:,:,secondegvali),2));
                    M(:,2) = squeeze(mean(incorrect_states(1,:,:,secondegvali),2));
                    S(:,1) = squeeze(std(correct_states(1,:,:,secondegvali),[],2));
                    S(:,2) = squeeze(std(incorrect_states(1,:,:,secondegvali),[],2));
                else
                    M=[];S=[];
                    v=[];
                    for ai=1:length(spon_states)
                        v(:,ai) = squeeze(correct_states{ai}(1,:,secondegvali));
                    end
                    M(:,1)=mean(v,2);
                    S(:,1)=std(v,[],2);
                    v=[];
                    for ai=1:length(spon_states)
                        v(:,ai) = squeeze(incorrect_states{ai}(1,:,secondegvali));
                    end
                    M(:,2)=mean(v,2);
                    S(:,2)=std(v,[],2);
                end
                subplot(3,1,2);
                plot_3_bars(M(:,1),S(:,1),statenames);ylim([0 .2]);title('c');
                subplot(3,1,3);plot_3_bars(M(:,2),S(:,2),statenames);ylim([0 .2]);title('i');
                mysave(gcf, fullfile(outputfiggolder,['second_eigval_'  isweigtedstr{isweigted}  '_bars_' signals_names{sig_i} '_th' num2str(th)]));
            end
            close all;
        end
        
    end
end
end
% spont
%     diffmask = scores_to_heatmap_allen(mean(spon_states_weighted.high_pup_l.(cent_features{ni})-spon_states_weighted.low_pup_q.(cent_features{ni}),2));
%     figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%     plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%         [],  -L, L, 1,colormap(redblue))
%     title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
%     % mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_weighted_heatmap_allen']));
%     % trial
%     for state_i=1:length(statenames)
%         diffmask = scores_to_heatmap_allen(mean(trial_states_weighted.(statenames{state_i}).(cent_features{ni}).correct-trial_states_weighted.(statenames{state_i}).(cent_features{ni}).incorrect,2));
%         figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%         plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%             [],  -L, L, 1,colormap(redblue))
%         title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
%
%         % mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_weighted_heatmap_allen']));
%     end
%     % not weighted
%     %     graph_overlay_allen_paired('',fullfile('', 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%     %         spon_states_notweighted.low_pup_q.(cent_features{ni}),strcat('','/spon_run_pupillow/'),cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spont)'],parcels_names,length(animals));
%     % spont
%     diffmask = scores_to_heatmap_allen(mean(spon_states_notweighted.high_pup_l.(cent_features{ni})-spon_states_notweighted.low_pup_q.(cent_features{ni}),2));
%     figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%     plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%         [],  -L, L, 1,colormap(redblue))
%     title(['run vs pupil low ' cent_features{ni} ' Centrality (spont)']);
%
%
%     % mysave(gcf, fullfile(outputfiggolder,['spont_state3_minus_state1_',cent_features{ni},'_notweighted_heatmap_allen']));
%     % trial not
%     for state_i=1:length(statenames)
%         diffmask = scores_to_heatmap_allen(mean(trial_states_notweighted.(statenames{state_i}).(cent_features{ni}).correct-trial_states_notweighted.(statenames{state_i}).(cent_features{ni}).incorrect,2));
%         figure;L=max(abs([min(diffmask(:)), max(diffmask(:))]));
%         plot_vals_heatmap(diffmask, 'Difference in Node Centrality',...
%             [],  -L, L, 1,colormap(redblue))
%         title(['Correct minus incorrect on ' statenames{state_i} ' ' cent_features{ni} ' Centrality (trial)']);
%         % mysave(gcf, fullfile(outputfiggolder,['trial_state' num2str(state_i) '_',cent_features{ni},'_notweighted_heatmap_allen']));
%     end


% end

% end




function plot_correct_incorrect_per_3parcels(M, S, parcels_names, statenames,spatialindex)

CondColors = get_3states_colors;
for s=1:length(statenames)
    statenames{s}(statenames{s}=='_') = ' ';
end
figure;
set(gcf,'renderer','painters');
for parcel_i = 1:length(spatialindex)
    subplot(length(spatialindex),1,(parcel_i));
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData(2:3)=0;
    h1(2).YData(2:3)=0;
    h1(1).FaceColor=CondColors(1,:);
    h1(2).FaceColor=CondColors(1,:);
    h1(2).FaceAlpha=0.4;
    h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');hold all;
    h1(1).YData([1 3])=0;
    h1(2).YData([1 3])=0;
    h1(1).FaceColor=CondColors(2,:);
    h1(2).FaceColor=CondColors(2,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    h1=barwitherr(squeeze(S(parcel_i,:,:))',squeeze(M(parcel_i,:,:))');
    h1(1).YData([1 2])=0;
    h1(2).YData([1 2])=0;
    h1(1).FaceColor=CondColors(3,:);
    h1(2).FaceColor=CondColors(3,:);
    h1(2).FaceAlpha=0.4;h1(1).LineStyle='none';
    h1(2).LineStyle='none';
    
    title(parcels_names{spatialindex(parcel_i)});
    
    set(gca,'xtick',1:23)
    set(gcf, 'Position',  [1,1, 700,1000]);
    set(gca,'xticklabel',statenames)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);%ylim([0 400]);
end
end


function makepsychbarplot(c50s,animals,statenames,name)
n = length(animals);
mean_mean_across_groups1=nanmean(c50s,1);
std_mean_across_groups1=nanstd(c50s,[],1)./sqrt(n-1);
figure;
CondColors = get_3states_colors;
subplot(1,1,1)
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, mean_mean_across_groups1(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:3)
set(gca,'xticklabel',statenames)
h = errorbar(1:3,mean_mean_across_groups1, std_mean_across_groups1,'LineStyle','none','LineWidth',0.5);title(name);
h.Color='k';
set(h, 'marker', 'none');
% mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\','figure_1',name), 'all');
end


%%
function eval_diff_map(animals, statenames)
addpath(genpath('../centrality_measures/'));
outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
for ai=1:length(animals)
    animal=animals{ai};
    
    for state_i = 1:length(statenames)
        load(fullfile(outputfolder,'network_analysis_corr',statenames{state_i}),'W_corr');
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
% function plot_centrality_res_gal(animals, outputfiggolder, statenames)
% [~, ~, finalindex] = get_allen_meta_parcels;
% cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
% for state_i = 1:length(statenames)
%     for cent_i = 1:length(cent_features)
%         spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = zeros(256,256,length(animals));
%     end
% end
% for i=1:length(animals)
%     animal=animals{i};
%     outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
%     for state_i = 1:length(statenames)
%         load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} 'LSSC']));
%         cent_features = fieldnames(cent_corr_notweighted);
%         for cent_i = 1:length(cent_features)
%             P = scores_to_heatmap_gal(cent_corr_notweighted.(cent_features{cent_i}), animal);
%             spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i})(:, :, i) = P;
%
%         end
%
%     end
% end
% mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
% for ni = 1:length(cent_features)-1
%     %difference maps, not weighted, for each centrality measure
%     braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
%     parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
%
%     diffmask = nanmean(spon_states_notweighted.high_pup_l.(cent_features{ni}),3)-...
%         nanmean(spon_states_notweighted.low_pup_q.(cent_features{ni}),3);
%     plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup loc minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
%     % mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_run_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
%     diffmask = nanmean(spon_states_notweighted.high_pup_q.(cent_features{ni}),3)-...
%         nanmean(spon_states_notweighted.low_pup_q.(cent_features{ni}),3);
%     plot_heatmap(diffmask, -0.04, 0.04, ['Difference in Centrality ' cent_features{ni} 'high pup q minus low pup q'], parcelsallen.parcells_new.indicators(:,:,finalindex));
%     % mysave(gcf, fullfile(outputfiggolder,'not_weighted','spon_pupilhigh_pupillow',strcat(cent_features{ni},'_heatmap_gal')), 'all');
% end
%
% end

% function plot_centrality_res(animals, outputfiggolder, statenames)
% [parcels_names] = get_allen_meta_parcels;
% cent_features = {'degree' 'closeness' 'betweenness' 'pagerank' 'eigenvector', 'participation', 'community'};
% for state_i = 1:length(statenames)
%     for cent_i = 1:length(cent_features)
%         spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = [];
%         %         spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = [];
%     end
% end
% for i=1:length(animals)
%     animal=animals{i};
%     outputfolder='X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\network_centrality';
%     for state_i = 1:length(statenames)
%         load(fullfile(outputfolder,[animal '_partial_corr_',statenames{state_i} 'Allen']));
%         cent_features = fieldnames(cent_corr_notweighted);
%         for cent_i = 1:length(cent_features)
%             spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
%                 spon_states_notweighted.(statenames{state_i}).(cent_features{cent_i}), ...
%                 cent_corr_notweighted.(cent_features{cent_i}));
%
%             %             spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}) = cat(2,...
%             %                 spon_states_weighted.(statenames{state_i}).(cent_features{cent_i}), ...
%             %                 cent_corr_weighted.(cent_features{cent_i}));
%         end
%
%     end
% end
% % mkNewDir(fullfile(outputfiggolder, 'weighted'))
% mkNewDir(fullfile(outputfiggolder, 'not_weighted'))
% legstr = {'Low Q', 'High Q', 'Loc'};
% for ni = 1:length(cent_features)-1
%     graph_overlay_allen_2conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         'spon_threestates',cent_features{ni},['2 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr([1 3]));
%
%
%
%     %     graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'weighted'), spon_states_weighted.low_pup_q.(cent_features{ni}),...
%     %         spon_states_weighted.high_pup_q.(cent_features{ni}), spon_states_weighted.high_pup_l.(cent_features{ni}),...
%     %         'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
%     %
%     graph_overlay_allen_3conditions([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         'spon_threestates',cent_features{ni},['3 states ' cent_features{ni} ' Centrality (spon)'], parcels_names,length(animals), legstr);
%
%     %difference maps, not weighted, for each centrality measure
%     braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
%     parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
%
%     %high vs low pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_q.(cent_features{ni}),...
%         spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         'spon_pupilhigh_pupillow',cent_features{ni},['pupil high vs low ' cent_features{ni} ' Centrality (spon)']);
%
%     %locomotion vs low pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         spon_states_notweighted.low_pup_q.(cent_features{ni}),'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.low_pup_q.(cent_features{ni}),...
%         'spon_run_pupillow',cent_features{ni},['run vs pupil low ' cent_features{ni} ' Centrality (spon)']);
%
%     %locomotion vs high pup
%     graph_overlay_allen_paired([],fullfile(outputfiggolder, 'not_weighted'), spon_states_notweighted.high_pup_l.(cent_features{ni}),...
%         spon_states_notweighted.high_pup_q.(cent_features{ni}),'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)'],parcels_names,length(animals));
%
%     graph_heatmap([],fullfile(outputfiggolder, 'not_weighted'),braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%         spon_states_notweighted.high_pup_l.(cent_features{ni}),spon_states_notweighted.high_pup_q.(cent_features{ni}),...
%         'spon_run_pupilhigh',cent_features{ni},['run vs pupil high ' cent_features{ni} ' Centrality (spon)']);
%
% end
%
% end
