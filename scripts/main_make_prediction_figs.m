% Figure 4
function main_make_prediction_figs
animals={'xw','xx','xz','xt','xu' 'xs'};
stateslabels = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
outputfiggolder = 'X:\Hadas\Meso-imaging\lan\meso_results\figs\';
%% Fig 3 - prediction_behavior
plot_prediction_acc(animals, stateslabels, fullfile(outputfiggolder, 'prediction_behavior'));
end



function plot_prediction_acc(animals, statenames, outputfiggolder)
mkNewDir(outputfiggolder);
 acc_per_parcel.Allen = nan(23, length(statenames), length(animals));
 acc_per_parcel.LSSC = [];
 acc_per_parcel.Grid4 = [];
 
signalsnames = {'Allen', 'LSSC', 'Grid4'};
for sig_i = 1:length(signalsnames)
    signal = signalsnames{sig_i};
    acc_all_parcels.(signal) = nan(length(statenames), length(animals));
    acc_heat.(signal) = nan(256,256,length(statenames), length(animals));
   
    for ai = 1:length(animals)
        for state_i = 1:length(statenames)
            predictionfile = fullfile('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory\behavior_prediction\',...
                [animals{ai} '_' statenames{state_i}, '_' signal '.mat']);
            if exist(predictionfile, 'file')                
                load(predictionfile,...
                    'accuracy_vec','accuracy_mat');
                acc_all_parcels.(signal)(state_i, ai) = mean(accuracy_vec);
                
                switch signal
                    case 'Allen'
                        acc_per_parcel.Allen(:, state_i, ai) = mean(accuracy_mat);
                        acc_heat.Allen(:,:,state_i, ai) = scores_to_heatmap_allen(mean(accuracy_mat)', false);
                        
                    case 'LSSC'
                        acc_per_parcel.LSSC = [];
                        acc_heat.LSSC(:,:,state_i, ai) = scores_to_heatmap_gal(mean(accuracy_mat)', animals{ai}, false);
                        
                    case 'Grid4'
                        acc_per_parcel.Grid4 = [];
                        acc_heat.Grid4(:,:,state_i, ai) = scores_to_heatmap_grid(mean(accuracy_mat)', animals{ai}, 4);
                end
            end
        end
    end
end
strstatenames=cell(3,1);
for k=1:length(statenames)
    strstatenames{k} =  statenames{k};
    strstatenames{k}(statenames{k}=='_') = ' ';
end
parcels_names.Allen = get_allen_meta_parcels;
n = length(animals);
% overall acc bars
for si=1:length(signalsnames)
    mean_mean_across_groups(:,si)=nanmean(acc_all_parcels.(signalsnames{si}),2);
    sem_mean_across_groups(:,si)=nanstd(acc_all_parcels.(signalsnames{si}),[],2)/sqrt(n-1);
end
figure;
plot_bars_3colors_per_parcel(mean_mean_across_groups',sem_mean_across_groups', strstatenames, signalsnames);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_all_parcels'));
% per parcel - Allen
M = nanmean(acc_per_parcel.Allen,3);
S = nanstd(acc_per_parcel.Allen,[],3)/sqrt(n-1);
plot_bars_3colors_per_parcel(M,S, strstatenames, parcels_names.Allen);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_per_parcel_Allen'));
% accuracy heat maps per state
figure;l=1;
for si=1:length(signalsnames)
    for k=1:length(statenames)
        subplot(3,3,l);
        %allen_parcels.indicators(:,:,finalindex.Allen)
        plot_vals_heatmap(nanmean(acc_heat.(signalsnames{si})(:,:,k,:),4), 'Prediction Accracy', [], 0.5, .7, true, colormap(jet));
        title([signalsnames{si} ' ' strstatenames{k}]);
        l=l+1;
    end
end
set(gcf,'Position',[680   239   864   739]);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_heatmaps'));

k=3;
% accuracy heat maps per animal on state 3
figure;l=1;
for si=1:length(signalsnames)
    for ai=1:length(animals)
        subplot(3,6,l);
        %allen_parcels.indicators(:,:,finalindex.Allen)
        plot_vals_heatmap(acc_heat.(signalsnames{si})(:,:,k,ai), 'Prediction Accracy', [], 0.5, .7, true, colormap(jet));
        title(['Animal #' num2str(ai)]);colormap jet;
        l=l+1;
    end
end
set(gcf,'Position',[116 301 1744 563]);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_animal_highpup_loc_heatmaps'));
% figure;
% for si=1:length(signals)
%     subplot(1,3,si);
%     plot_vals_heatmap(mean(acc_heat.(signals{si})(:,:,3,:)-acc_heat.(signals{si})(:,:,1,:),4), 'Prediction Accracy', [], -0.1, .1, true, colormap(redblue));
% end

% graph_overlay_allen_paired([],outputfiggolder,M(:,3),M(:,1),'','run_pupillow_accuracy_SVM_diff','run - pupil low (svm acc)',parcels_names,n)
% graph_heatmap([],outputfiggolder,braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%     M(:,3),M(:,1),'','run_pupillow_accuracy_SVM_heatmap','run - pupil low (svm acc)');
%
% graph_overlay_allen_paired([],'not_weighted',M(:,3),M(:,2),'','run_pupilhigh_accuracy_SVM_diff','run - pupil high (svm acc)',parcels_names,n)
% graph_heatmap([],outputfiggolder,braininfo.brain_mask,parcelsallen.parcells_new.indicators,...
%     M(:,3),M(:,2),'','run_pupilhigh_accuracy_SVM_heatmap','run - pupil high (svm acc)');
end



