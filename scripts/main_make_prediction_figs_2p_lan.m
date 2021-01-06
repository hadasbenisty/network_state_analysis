% Figure 4
function main_make_prediction_figs_2p_lan
T = readtable('X:\Lan\FMB208 conditioning imaging data\animal_log_imaging_lt.xlsx');
animals = T.AnimalID;

poptype = T.TargetPop;
stateslabels = {'low_pup_q', 'high_pup_q', 'high_pup_l'};
outputfiggolder = 'X:\Hadas\Meso-imaging\lan\meso_results\figs2p\';
mkNewDir(outputfiggolder);
%% Fig 3 - prediction_behavior

plot_prediction_acc(animals, poptype, stateslabels, fullfile(outputfiggolder, 'prediction_behavior'));

end



function plot_prediction_acc(animals, poptype, statenames, outputfiggolder)
mkNewDir(outputfiggolder);

acc_all_parcels = nan(length(statenames), length(animals));
acc_per_cell = cell(length(animals),1);

for ai = 1:length(animals)
    for state_i = 1:length(statenames)
        predictionfile = fullfile('X:\Hadas\Meso-imaging\lan\meso_results\ProcessingDirectory2p\behavior_prediction\',...
            [animals{ai} '_' statenames{state_i} '.mat']);
        if exist(predictionfile, 'file')
            load(predictionfile,...
                'accuracy_vec','accuracy_mat');
            acc_all_parcels(state_i, ai) = mean(accuracy_vec);
            
            
            acc_per_cell{ai}(:, state_i) = mean(accuracy_mat);
             acc_per_cell_std{ai}(:, state_i) = std(accuracy_mat);
                
        end
    end

end
strstatenames=cell(3,1);
for k=1:length(statenames)
    strstatenames{k} =  statenames{k};
    strstatenames{k}(statenames{k}=='_') = ' ';
end
populations = unique(poptype);
% overall acc bars per population

for pi=1:length(populations)
    curpop = populations{pi};
    curraninmals = (strcmp(poptype, curpop));
n = sum(~isnan(acc_all_parcels(:,curraninmals)),2);

mean_mean_across_groups(:,pi)=nanmean(acc_all_parcels(:,curraninmals),2);
sem_mean_across_groups(:,pi)=nanstd(acc_all_parcels(:,curraninmals),[],2)./sqrt(n-1);
end
figure;
plot_bars_3colors(mean_mean_across_groups',sem_mean_across_groups', strstatenames,populations);
mysave(gcf,fullfile(outputfiggolder, 'behavior_prediction_by_state_all_cells'));
% per nrn
binsN = 64;
bins = linspace(0, 1, binsN);
data4hist = cell(length(statenames),length(populations));
for ai=1:length(animals)
    curpop = find(strcmp(populations,poptype{ai}));
    if ~isempty(acc_per_cell{ai})
    for si=1:size(acc_per_cell{ai},2)
        if ~all( acc_per_cell{ai}(:, si)==0)
    data4hist{si, curpop} = cat(1, data4hist{si, curpop}, acc_per_cell{ai}(:, si));
        end
    end
    end
end
figure;
Colors = get_3states_colors;
for pi = 1:length(populations)
        subplot(1,3,pi);
        for state_i = 1:length(statenames)
        
        histstates(:, state_i, pi) = hist(data4hist{state_i, pi}, bins);
       histstates(:, state_i, pi)=histstates(:, state_i, pi)/sum(histstates(:, state_i, pi));
        b=bar(bins, histstates(:, state_i, pi));
        b.FaceColor = Colors(state_i,:);
        hold all;
        end
    title(populations{pi});
    xlabel('Accuracy');
end
set(gcf,'Position', [1 41 1920 963]);
mysave(gcf,fullfile(outputfiggolder, 'acc_hist_nrns_populations'));

figure;l=1;
for pi = 1:length(populations)
        for state_i = 1:length(statenames)
            for state_j = 1:length(statenames)
     [~,ksres_bypop(state_i, state_j, pi)] = kstest2(data4hist{state_i, pi},data4hist{state_j, pi});
            end
        end
        subplot(1,3,l);
        imagesc(ksres_bypop(:,:,pi)<0.05);colorbar;l=l+1;
        title(populations{pi});
        set(gca,'XTick',1:length(statenames))
        set(gca,'XTickLabel',statenames);
        set(gca,'YTick',1:length(statenames))
        set(gca,'YTickLabel',statenames);
        
end


% for state_i = 1:length(statenames)
% for pi = 1:length(populations)
%      for pj = 1:length(populations)
%    
%             for state_j = 1:length(statenames)
%      [~,ksres_bystate(pi,pj,state_i)] = kstest2(data4hist{state_i, pi},data4hist{state_i, pj});
%             end
%         end
% end
% subplot(2,3,l);
%         imagesc(ksres_bystate(:,:,state_i)<0.05);colorbar;l=l+1;
%         title(statenames{state_i});
%         set(gca,'XTick',1:length(populations))
%         set(gca,'XTickLabel',populations);
%         set(gca,'YTick',1:length(populations))
%         set(gca,'YTickLabel',populations);
% end
suptitle(sprintf('KS test comparing histograms of prediction accuracies\nYellow = a significant difference'))
mysave(gcf,fullfile(outputfiggolder, 'kstest_nrns_populations'));

for ai = 1:length(animals)
    if ~isempty(acc_per_cell{ai})
   plot_bars_3colors_per_parcel(acc_per_cell{ai}, acc_per_cell_std{ai}, statenames(1:size(acc_per_cell{ai},2)), 1:size(acc_per_cell_std{ai},1));
    suptitle(['Single Cell Prediction ' animals{ai} ' ' poptype{ai}]);
    mysave(gcf,fullfile(outputfiggolder, ['single_cell_acc_bystate_' animals{ai} '_' poptype{ai}]));

    end
end
% histData = nan(binsN, length(statenames), length(animals));
% 
% for ai=1:length(animals)
%     if ~isempty(acc_per_cell{ai})
%         for si=1:length(statenames)
%             histData(:, si, ai) = hist(acc_per_cell{ai}(:,si), bins);
%         end
%     end
% end
end



