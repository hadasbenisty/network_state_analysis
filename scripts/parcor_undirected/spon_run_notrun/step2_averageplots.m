function step2_averageplots(inputpth, outputpth)
%% average plots for partial corr for spon run not run
%clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
spon_run_cat_degree=[];spon_notrun_cat_degree=[];
spon_run_cat_closeness=[];spon_notrun_cat_closeness=[];
spon_run_cat_betweenness=[];spon_notrun_cat_betweenness=[];
spon_run_cat_pagerank=[];spon_notrun_cat_pagerank=[];
spon_run_cat_eigenvector=[];spon_notrun_cat_eigenvector=[];
condition1='run';condition2='notrun';
parcels_region_labels_binotrunral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_binotrunral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    spon_run=load(strcat(inputpth,animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    spon_notrun=load(strcat(inputpth,animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    spon_run_cat_degree=cat(2, spon_run_cat_degree,spon_run.cent_corr.degree); 
    spon_run_cat_closeness=cat(2, spon_run_cat_closeness,spon_run.cent_corr.closeness); 
    spon_run_cat_betweenness=cat(2, spon_run_cat_betweenness,spon_run.cent_corr.betweenness);    
    spon_run_cat_pagerank=cat(2, spon_run_cat_pagerank,spon_run.cent_corr.pagerank); 
    spon_run_cat_eigenvector=cat(2, spon_run_cat_eigenvector,spon_run.cent_corr.eigenvector); 
     subplot(2,length(animals),counter)
    p=plot(spon_run.G_corr);title(strcat(condition1,animal));counter=counter+1;
    p.NodeCData=parcels_region_labels;
    p.NodeLabel = {};
    subplot(2,length(animals),counter)
    n=plot(spon_notrun.G_corr);title(strcat(condition2,animal));
    n.NodeCData=parcels_region_labels;counter=counter+1;
    n.NodeLabel = {};
    %mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_spon_network\',animal,'\','run_notrun_graph'),'all')
    
    spon_notrun_cat_degree=cat(2, spon_notrun_cat_degree,spon_notrun.cent_corr.degree); 
    spon_notrun_cat_closeness=cat(2, spon_notrun_cat_closeness,spon_notrun.cent_corr.closeness); 
    spon_notrun_cat_betweenness=cat(2, spon_notrun_cat_betweenness,spon_notrun.cent_corr.betweenness);       
    spon_notrun_cat_pagerank=cat(2, spon_notrun_cat_pagerank,spon_notrun.cent_corr.pagerank); 
    spon_notrun_cat_eigenvector=cat(2, spon_notrun_cat_eigenvector,spon_notrun.cent_corr.eigenvector); 
    clearvars spon_run spon_notrun
end
%get parcel label
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
parcels_names=textt(finalindex,1);

graph_overlay_allen_paired(outputpth, spon_run_cat_degree,spon_notrun_cat_degree,'spon_run_notrun','degree_centrality','run vs notrun degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth,spon_run_cat_closeness,spon_notrun_cat_closeness,'spon_run_notrun','closeness_centrality','run vs notrun closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth,spon_run_cat_pagerank,spon_notrun_cat_pagerank,'spon_run_notrun','pagerank_centrality','run vs notrun pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth,spon_run_cat_eigenvector,spon_notrun_cat_eigenvector,'spon_run_notrun','eigenvector_centrality','run vs notrun eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth,spon_run_cat_betweenness,spon_notrun_cat_betweenness,'spon_run_notrun','betweenness_centrality','run vs notrun betweenness Centrality (spon)',parcels_names,length(animals));

graph_overlay_allen_notpaired(outputpth, spon_run_cat_degree,spon_notrun_cat_degree,'spon_run_notrun','degree_centrality','run vs notrun degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth,spon_run_cat_closeness,spon_notrun_cat_closeness,'spon_run_notrun','closeness_centrality','run vs notrun closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth,spon_run_cat_pagerank,spon_notrun_cat_pagerank,'spon_run_notrun','pagerank_centrality','run vs notrun pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth,spon_run_cat_eigenvector,spon_notrun_cat_eigenvector,'spon_run_notrun','eigenvector_centrality','run vs notrun eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth,spon_run_cat_betweenness,spon_notrun_cat_betweenness,'spon_run_notrun','betweenness_centrality','run vs notrun betweenness Centrality (spon)',parcels_names,length(animals));

