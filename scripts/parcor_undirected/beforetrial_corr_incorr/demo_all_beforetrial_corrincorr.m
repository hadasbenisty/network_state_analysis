function demo_all_beforetrial_corrincorr(outputpth, animal, isw, saveplots)
disp(animal)
if nargin == 3
    saveplots = false;
end
rehash
rehash path
addpath(genpath('X:\Lav/network_state_analysis\functions'))
cd('X:\Lav\ProcessingDirectory')
inputfolder=fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',animal,'\');
outputfolder=fullfile(outputpth, animal,'\');
mkNewDir(outputfolder);
% load the time traces file (notice the animal and day)
rundata=load(strcat(inputfolder, strcat(animal,'state_trials_networkanalysis.mat')));
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_region_labels_bilateral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bilateral(finalindex);
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
parcels_names=textt(finalindex,1);

for j=1:10
    switch j
        case 1
            data = rundata.Correct(finalindex,all(~isnan(rundata.Correct)));
            condition='trial_correct';
        case 2
            data = rundata.Incorrect(finalindex,all(~isnan(rundata.Incorrect)));
            condition='trial_incorrect';
        case 3
            data = rundata.Correct_HighPupZscored(finalindex,all(~isnan(rundata.Correct_HighPupZscored)));
            condition='trial_highpup_correct';
        case 4
            data = rundata.Incorrect_HighPupZscored(finalindex,all(~isnan(rundata.Incorrect_HighPupZscored)));
            condition='trial_highpup_incorrect';
        case 5
            data = rundata.Correct_LowPupZscored(finalindex,all(~isnan(rundata.Correct_LowPupZscored)));
            condition='trial_lowpup_correct';
        case 6
            data = rundata.Incorrect_LowPupZscored(finalindex,all(~isnan(rundata.Incorrect_LowPupZscored)));
            condition='trial_lowpup_incorrect';
        case 7
            data = rundata.Correct_NotRunZscored(finalindex,all(~isnan(rundata.Correct_NotRunZscored)));
            condition='trial_notrun_correct';   % for nan - columns (so all parcels are aligned).
        case 8
            data = rundata.Incorrect_NotRunZscored(finalindex,all(~isnan(rundata.Incorrect_NotRunZscored)));
            condition='trial_notrun_incorrect';   % for nan - columns (so all parcels are aligned).
        case 9
            data = rundata.Incorrect_RunZscored(finalindex,all(~isnan(rundata.Incorrect_RunZscored)));
            condition='trial_run_incorrect';   % for nan - columns (so all parcels are aligned).
        case 10
            data = rundata.Correct_RunZscored(finalindex,all(~isnan(rundata.Correct_RunZscored)));
            condition='trial_run_correct';   % for nan - columns (so all parcels are aligned).
    end
    disp(condition)
    if isempty(find(isnan(data)==1))
        %% measure weights
        W_corr = measure_weights_partial(data, 'corr');
        disp('corweights done')
        %         [W_lasso, well_modeled_nodes] = measure_weights(data, 'relative_modeling_contribution');
        %         disp('lassoweights done')
        
        %% Graph Analysis
        % get threshold for connected nodes
        
        [indic_corr, cent_corr, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names, isw);
        %         [indic_lasso, cent_lasso, G_lasso, names_lasso] = graph_analysis_afterclust(W_lasso, parcels_names);
        disp('graph analysis done')
        
        save(strcat(outputfolder,'network_analysis_corr',condition),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
        %         save(strcat(outputfolder,'network_analysis_lasso',condition),'W_lasso','well_modeled_nodes','indic_lasso','cent_lasso','G_lasso','names_lasso');
        disp('graph analysis saved')
        
        %% Visualization
        if saveplots
            % plot graphs
            plot_graph(G_corr, cent_corr, names_corr, parcels_region_labels);
            suptitle('Correlation');
            %set(gcf, 'Position',  [150,150, 1000,500]);
            mysave(gcf, strcat(outputfolder,'correlation_graph',condition), 'fig');
            %
            %        plot_graph(G_lasso, cent_lasso, names_lasso, parcels_region_labels);
            %        suptitle('LASSO');
            %        set(gcf, 'Position',  [150,150, 1000,500]);
            %        mysave(gcf, strcat(outputfolder,'lasso_graph',condition), 'fig');
            %        disp('graphs plotted and saved')
            
            % plot centrality by nodes
            figure;
            for k = 1:length(names_corr)
                subplot(2,4,k);
                bar(cent_corr.(names_corr{k}));
                set(gca, 'XTickLabel', parcels_names);
                set(gca,'XTickLabelRotation',45);
                %set(gcf, 'Position',  [150,150, 1000,500]);
                title(names_corr{k});
            end
            suptitle('Correlation');
            mysave(gcf, strcat(outputfolder,'correlation_centrality',condition), 'fig');
            disp('cntrality plotted and saved')
            
            %
            %        figure;
            %        for k = 1:length(names_lasso)
            %            subplot(2,4,k);
            %            bar(cent_corr.(names_lasso{k}));
            %            set(gca, 'XTickLabel', parcels_names);
            %            set(gca,'XTickLabelRotation',45);
            %            set(gcf, 'Position',  [150,150, 1000,500]);
            %            title(names_lasso{k});
            %        end
            %        suptitle('LASSO');
            %        mysave(gcf, strcat(outputfolder,'lasso_centrality',condition), 'fig');
            %
            % plot centrality by node population
            figure;
            bar(sum(indic_corr,2)/size(indic_corr, 2));set(gca,'XTickLabel',names_corr);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
            title('Correlation');
            mysave(gcf, strcat(outputfolder,'centrality_by_node_pop',condition), 'fig');
        end
        disp(strcat('done',condition))
    else
        disp('nans in dataset')
    end
    clearvars data W_corr indic_corr cent_corr G_corr names_corr %W_lasso well_modeled_nodes indic_lasso cent_lasso G_lasso names_lasso
end

disp('done')

end

