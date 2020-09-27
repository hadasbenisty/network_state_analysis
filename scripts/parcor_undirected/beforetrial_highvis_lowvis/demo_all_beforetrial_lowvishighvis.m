function demo_all_beforetrial_lowvishighvis(animal)
rehash
rehash path
addpath(genpath('X:\Lav/network_state_analysis\functions'))
cd('X:\Lav\ProcessingDirectory')
inputfolder=fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude',animal,'\');
outputfolder=fullfile('X:\Lav\ProcessingDirectory\',animal,'\');
mkdir(animal);
% load the time traces file (notice the animal and day)
beforetrialdata=load(strcat(inputfolder, strcat(animal,'before_trials_networkanalysis.mat')));
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_region_labels_bilateral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bilateral(finalindex);
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
parcels_names=textt(finalindex,1);

for j=1:10
    if j==1
        data = beforetrialdata.HighVis(finalindex,all(~isnan(beforetrialdata.HighVis)));condition='trial_HighVis';
    elseif j==2
        data = beforetrialdata.LowVis(finalindex,all(~isnan(beforetrialdata.LowVis)));condition='trial_LowVis';
    elseif j==3
        data = beforetrialdata.HighVis_HighPupZscored(finalindex,all(~isnan(beforetrialdata.HighVis_HighPupZscored)));condition='trial_highpup_HighVis';
        %data = rundata.HighPupZscored(finalindex,all(~isnan(rundata.HighPupZscored)));condition='trial_highpup';   % for nan - columns (so all parcels are aligned).
    elseif j==4
        data = beforetrialdata.LowVis_HighPupZscored(finalindex,all(~isnan(beforetrialdata.LowVis_HighPupZscored)));condition='trial_highpup_LowVis';
        %data = rundata.LowPupZscored(finalindex,all(~isnan(rundata.LowPupZscored)));condition='trial_lowpup';      % for nan - columns (so all parcels are aligned).
    elseif j==5
        data = beforetrialdata.HighVis_LowPupZscored(finalindex,all(~isnan(beforetrialdata.HighVis_LowPupZscored)));condition='trial_lowpup_HighVis';
        %data = rundata.HighVis_RunZscored(finalindex,all(~isnan(rundata.HighVis_RunZscored)));condition='trial_run_HighVis';   % for nan - columns (so all parcels are aligned).
    elseif j==6
        data = beforetrialdata.LowVis_LowPupZscored(finalindex,all(~isnan(beforetrialdata.LowVis_LowPupZscored)));condition='trial_lowpup_LowVis';
        %data = rundata.LowVis_RunZscored(finalindex,all(~isnan(rundata.LowVis_RunZscored)));condition='trial_run_LowVis';   % for nan - columns (so all parcels are aligned).
    elseif j==7
        data = beforetrialdata.HighVis_NotRunZscored(finalindex,all(~isnan(beforetrialdata.HighVis_NotRunZscored)));condition='trial_notrun_HighVis';   % for nan - columns (so all parcels are aligned).
    elseif j==8
        data = beforetrialdata.LowVis_NotRunZscored(finalindex,all(~isnan(beforetrialdata.LowVis_NotRunZscored)));condition='trial_notrun_LowVis';   % for nan - columns (so all parcels are aligned).
    elseif j==9
        data = beforetrialdata.LowVis_RunZscored(finalindex,all(~isnan(beforetrialdata.LowVis_RunZscored)));condition='trial_run_LowVis';   % for nan - columns (so all parcels are aligned).
    elseif j==10
        data = beforetrialdata.HighVis_RunZscored(finalindex,all(~isnan(beforetrialdata.HighVis_RunZscored)));condition='trial_run_HighVis';   % for nan - columns (so all parcels are aligned).
    end
    if isempty(find(isnan(data)==1))
        %% measure weights
        W_corr = measure_weights_partial(data, 'corr');
        disp('corweights done')
        
        %% Graph Analysis
        % get threshold for connected nodes
        [indic_corr, cent_corr, G_corr, names_corr] = graph_analysis_afterclust(W_corr, parcels_names);
        disp('graph analysis done')
        
        save(strcat(outputfolder,'network_analysis_corr',condition),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
        disp('graph analysis saved')
        
        %% Visualization
        % plot graphs
        plot_graph(G_corr, cent_corr, names_corr, parcels_region_labels);
        suptitle('Correlation');
        set(gcf, 'Position',  [150,150, 1000,500]);
        mysave(gcf, strcat(outputfolder,'correlation_graph',condition), 'fig');
        
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
        mysave(gcf, strcat(outputfolder,'correlation_centrality',condition), 'fig');
        disp('cntrality plotted and saved')
        % plot centrality by node population
        figure;
        bar(sum(indic_corr,2)/size(indic_corr, 2));set(gca,'XTickLabel',names_corr);set(gca,'XTickLabelRotation',45);set(gcf, 'Position',  [150,150, 1000,500]);
        title('Correlation');
        mysave(gcf, strcat(outputfolder,'centrality_by_node_pop',condition), 'fig');
        disp(strcat('done',condition))
    else
        disp('nans in dataset')
    end
    clearvars data W_corr indic_corr cent_corr G_corr names_corr
end

disp('done')

end

