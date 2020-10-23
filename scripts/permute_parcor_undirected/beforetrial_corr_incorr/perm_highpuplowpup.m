function perm_beforetrial_corincorr(animal)
rehash
rehash path
addpath(genpath('X:\Lav/network_state_analysis\functions'))
cd('X:\Lav\ProcessingDirectory')
inputfolder=fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',animal,'\');
outputfolder=fullfile('X:\Lav\ProcessingDirectory\',animal,'\permute\');
mkNewDir(animal);
% load the time traces file (notice the animal and day)
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_region_labels_bilateral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bilateral(finalindex);
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
parcels_names=textt(finalindex,1);
for kr=1:100
    % load the time traces file (notice the animal and day)
    rundata=load(strcat(inputfolder, strcat(num2str(kr),animal,'state_trials_networkanalysis.mat')));
    for j=1:10
        if j==1
            data = rundata.Correct(finalindex,all(~isnan(rundata.Correct)));condition='trial_correct';
        elseif j==2
            data = rundata.Incorrect(finalindex,all(~isnan(rundata.Incorrect)));condition='trial_incorrect';
        elseif j==3
            data = rundata.Correct_HighPupZscored(finalindex,all(~isnan(rundata.Correct_HighPupZscored)));condition='trial_highpup_correct';
        elseif j==4
            data = rundata.Incorrect_HighPupZscored(finalindex,all(~isnan(rundata.Incorrect_HighPupZscored)));condition='trial_highpup_incorrect';
        elseif j==5
            data = rundata.Correct_LowPupZscored(finalindex,all(~isnan(rundata.Correct_LowPupZscored)));condition='trial_lowpup_correct';
        elseif j==6
            data = rundata.Incorrect_LowPupZscored(finalindex,all(~isnan(rundata.Incorrect_LowPupZscored)));condition='trial_lowpup_incorrect';
        elseif j==7
            data = rundata.Correct_NotRunZscored(finalindex,all(~isnan(rundata.Correct_NotRunZscored)));condition='trial_notrun_correct';   % for nan - columns (so all parcels are aligned).
        elseif j==8
            data = rundata.Incorrect_NotRunZscored(finalindex,all(~isnan(rundata.Incorrect_NotRunZscored)));condition='trial_notrun_incorrect';   % for nan - columns (so all parcels are aligned).
        elseif j==9
            data = rundata.Incorrect_RunZscored(finalindex,all(~isnan(rundata.Incorrect_RunZscored)));condition='trial_run_incorrect';   % for nan - columns (so all parcels are aligned).
        elseif j==10
            data = rundata.Correct_RunZscored(finalindex,all(~isnan(rundata.Correct_RunZscored)));condition='trial_run_correct';   % for nan - columns (so all parcels are aligned).
        end
        if isempty(find(isnan(data)==1))
            %% measure weights
            W_corr = measure_weights_partial(data, 'corr');
            disp('corweights done')
            %[W_lasso, well_modeled_nodes] = measure_weights(data, 'relative_modeling_contribution');
            %disp('lassoweights done')
            %% Graph Analysis
            % get threshold for connected nodes
            [indic_corr, cent_corr, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names);
            %[indic_lasso, cent_lasso, G_lasso, names_lasso] = graph_analysis_afterclust(W_lasso, parcels_names);
            
            disp('graph analysis done')
            save(strcat(outputfolder,num2str(kr),'network_analysis_corr',condition),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
            %save(strcat(outputfolder,num2str(kr),'network_analysis_lasso',condition),'W_lasso','indic_lasso','cent_lasso','G_lasso','names_lasso');
            disp('saved')
            if kr==1|kr==50
                %% Visualization
                % plot graphs
                plot_graph(G_corr, cent_corr, names_corr, parcels_region_labels);
                suptitle('Correlation');
                set(gcf, 'Position',  [150,150, 1000,500]);
                mysave(gcf, strcat(outputfolder,num2str(kr),'correlation_graph',condition), 'fig');
                
                %                 plot_graph(G_lasso, cent_lasso, names_lasso, parcels_region_labels);
                %                 suptitle('LASSO');
                %                 set(gcf, 'Position',  [150,150, 1000,500]);
                %                 mysave(gcf, strcat(outputfolder,num2str(kr),'lasso_graph',condition), 'fig');
                %                 disp('graphs plotted and saved')
                
                % plot centrality by nodes
                figure;
                for k = 1:length(names_corr)
                    subplot(2,4,k);
                    bar(cent_corr.(names_corr{k}));
                    set(gca, 'XTickLabel', parcels_names);
                    set(gca,'XTickLabelRotation',45);
                    set(gcf, 'Position',  [150,150, 1000,500]);
                    title(names_corr{k});
                end
                suptitle('Correlation');
                mysave(gcf, strcat(outputfolder,num2str(kr),'correlation_centrality',condition), 'fig');
                disp('cntrality plotted and saved')
                
                %
                %                 figure;
                %                 for k = 1:length(names_lasso)
                %                     subplot(2,4,k);
                %                     bar(cent_corr.(names_lasso{k}));
                %                     set(gca, 'XTickLabel', parcels_names);
                %                     set(gca,'XTickLabelRotation',45);
                %                     set(gcf, 'Position',  [150,150, 1000,500]);
                %                     title(names_lasso{k});
                %                 end
                %                 suptitle('LASSO');
                %                 mysave(gcf, strcat(outputfolder,num2str(kr),'lasso_centrality',condition), 'fig');
                %
                % plot centrality by node population
                figure;
                subplot(2,1,1);
                bar(sum(indic_corr,2)/size(indic_corr, 2));set(gca,'XTickLabel',names_corr);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
                title('Correlation');
                %                 subplot(2,1,2);
                %                 bar(sum(indic_lasso(:, well_modeled_nodes),1)/length(well_modeled_nodes));set(gca,'XTickLabel',names_lasso);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
                %                 title('LASSO');
                %                 suptitle('Centrality Stats of Well Modeled Nodes');
                mysave(gcf, strcat(outputfolder,num2str(kr),'centrality_by_node_pop',condition), 'fig');
            else
            end
            
        else
            disp('nans in dataset')
        end
        clearvars data W_corr indic_corr cent_corr G_corr names_corr %W_lasso well_modeled_nodes
    end
    disp(strcat('done',num2str(kr)));
end
end

