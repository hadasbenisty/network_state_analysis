%% corr
clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
trial_notrun_highvis_cat_degree=[];trial_notrun_lowvis_cat_degree=[];
trial_notrun_highvis_cat_closeness=[];trial_notrun_lowvis_cat_closeness=[];
trial_notrun_highvis_cat_betweenness=[];trial_notrun_lowvis_cat_betweenness=[];
trial_notrun_highvis_cat_pagerank=[];trial_notrun_lowvis_cat_pagerank=[];
trial_notrun_highvis_cat_eigenvector=[];trial_notrun_lowvis_cat_eigenvector=[];
condition1='trial_notrun_highvis';condition2='trial_notrun_lowvis';
parcels_region_labels_bitrial_notrun_lowvisral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bitrial_notrun_lowvisral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    trial_notrun_highvis=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    trial_notrun_lowvis=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    trial_notrun_highvis_cat_degree=cat(2, trial_notrun_highvis_cat_degree,trial_notrun_highvis.cent_corr.degree); 
    trial_notrun_highvis_cat_closeness=cat(2, trial_notrun_highvis_cat_closeness,trial_notrun_highvis.cent_corr.closeness); 
    trial_notrun_highvis_cat_betweenness=cat(2, trial_notrun_highvis_cat_betweenness,trial_notrun_highvis.cent_corr.betweenness);    
    trial_notrun_highvis_cat_pagerank=cat(2, trial_notrun_highvis_cat_pagerank,trial_notrun_highvis.cent_corr.pagerank); 
    trial_notrun_highvis_cat_eigenvector=cat(2, trial_notrun_highvis_cat_eigenvector,trial_notrun_highvis.cent_corr.eigenvector); 
     subplot(2,length(animals),counter)
    p=plot(trial_notrun_highvis.G_corr);title(strcat(condition1,animal));counter=counter+1;
    p.NodeCData=parcels_region_labels;
    p.NodeLabel = {};
    subplot(2,length(animals),counter)
    n=plot(trial_notrun_lowvis.G_corr);title(strcat(condition2,animal));
    n.NodeCData=parcels_region_labels;counter=counter+1;
    n.NodeLabel = {};
    %mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_network\',animal,'\','trial_notrun_highvis_trial_notrun_lowvis_graph'),'all')
    
    trial_notrun_lowvis_cat_degree=cat(2, trial_notrun_lowvis_cat_degree,trial_notrun_lowvis.cent_corr.degree); 
    trial_notrun_lowvis_cat_closeness=cat(2, trial_notrun_lowvis_cat_closeness,trial_notrun_lowvis.cent_corr.closeness); 
    trial_notrun_lowvis_cat_betweenness=cat(2, trial_notrun_lowvis_cat_betweenness,trial_notrun_lowvis.cent_corr.betweenness);       
    trial_notrun_lowvis_cat_pagerank=cat(2, trial_notrun_lowvis_cat_pagerank,trial_notrun_lowvis.cent_corr.pagerank); 
    trial_notrun_lowvis_cat_eigenvector=cat(2, trial_notrun_lowvis_cat_eigenvector,trial_notrun_lowvis.cent_corr.eigenvector); 
    clearvars trial_notrun_highvis trial_notrun_lowvis
end
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_names=textt(finalindex,1);

graph_overlay_allen_paired(trial_notrun_highvis_cat_degree,trial_notrun_lowvis_cat_degree,'trial_notrun_highvis_trial_notrun_lowvis','degree_centrality','trial_notrun_highvis vs trial_notrun_lowvis degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(trial_notrun_highvis_cat_closeness,trial_notrun_lowvis_cat_closeness,'trial_notrun_highvis_trial_notrun_lowvis','closeness_centrality','trial_notrun_highvis vs trial_notrun_lowvis closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(trial_notrun_highvis_cat_pagerank,trial_notrun_lowvis_cat_pagerank,'trial_notrun_highvis_trial_notrun_lowvis','pagerank_centrality','trial_notrun_highvis vs trial_notrun_lowvis pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(trial_notrun_highvis_cat_eigenvector,trial_notrun_lowvis_cat_eigenvector,'trial_notrun_highvis_trial_notrun_lowvis','eigenvector_centrality','trial_notrun_highvis vs trial_notrun_lowvis eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(trial_notrun_highvis_cat_betweenness,trial_notrun_lowvis_cat_betweenness,'trial_notrun_highvis_trial_notrun_lowvis','betweenness_centrality','trial_notrun_highvis vs trial_notrun_lowvis betweenness Centrality (spon)',parcels_names,length(animals));

graph_overlay_allen_notpaired(trial_notrun_highvis_cat_eigenvector,trial_notrun_lowvis_cat_eigenvector,'trial_notrun_highvis_trial_notrun_lowvis','eigenvector_centrality','trial_notrun_highvis vs trial_notrun_lowvis eigenvector Centrality (spon)',parcels_names,length(animals));
