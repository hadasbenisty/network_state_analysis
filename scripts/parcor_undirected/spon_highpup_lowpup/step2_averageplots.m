function step2_averageplots(inputpth, outputpth)
%% corr
% clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
spon_pupilhigh_cat_degree=[];spon_pupillow_cat_degree=[];
spon_pupilhigh_cat_closeness=[];spon_pupillow_cat_closeness=[];
spon_pupilhigh_cat_betweenness=[];spon_pupillow_cat_betweenness=[];
spon_pupilhigh_cat_pagerank=[];spon_pupillow_cat_pagerank=[];
spon_pupilhigh_cat_eigenvector=[];spon_pupillow_cat_eigenvector=[];
condition1='pupilhigh';condition2='pupillow';
parcels_region_labels_bipupillowral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bipupillowral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    spon_pupilhigh=load(strcat(inputpth,animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    spon_pupillow=load(strcat(inputpth,animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    spon_pupilhigh_cat_degree=cat(2, spon_pupilhigh_cat_degree,spon_pupilhigh.cent_corr.degree); 
    spon_pupilhigh_cat_closeness=cat(2, spon_pupilhigh_cat_closeness,spon_pupilhigh.cent_corr.closeness); 
    spon_pupilhigh_cat_betweenness=cat(2, spon_pupilhigh_cat_betweenness,spon_pupilhigh.cent_corr.betweenness);    
    spon_pupilhigh_cat_pagerank=cat(2, spon_pupilhigh_cat_pagerank,spon_pupilhigh.cent_corr.pagerank); 
    spon_pupilhigh_cat_eigenvector=cat(2, spon_pupilhigh_cat_eigenvector,spon_pupilhigh.cent_corr.eigenvector); 
     subplot(2,length(animals),counter)
    p=plot(spon_pupilhigh.G_corr);title(strcat(condition1,animal));counter=counter+1;
    p.NodeCData=parcels_region_labels;
    p.NodeLabel = {};
    subplot(2,length(animals),counter)
    n=plot(spon_pupillow.G_corr);title(strcat(condition2,animal));
    n.NodeCData=parcels_region_labels;counter=counter+1;
    n.NodeLabel = {};
    %mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_spon_network\',animal,'\','pupilhigh_pupillow_graph'),'all')
    
    spon_pupillow_cat_degree=cat(2, spon_pupillow_cat_degree,spon_pupillow.cent_corr.degree); 
    spon_pupillow_cat_closeness=cat(2, spon_pupillow_cat_closeness,spon_pupillow.cent_corr.closeness); 
    spon_pupillow_cat_betweenness=cat(2, spon_pupillow_cat_betweenness,spon_pupillow.cent_corr.betweenness);       
    spon_pupillow_cat_pagerank=cat(2, spon_pupillow_cat_pagerank,spon_pupillow.cent_corr.pagerank); 
    spon_pupillow_cat_eigenvector=cat(2, spon_pupillow_cat_eigenvector,spon_pupillow.cent_corr.eigenvector); 
    clearvars spon_pupilhigh spon_pupillow
end
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_names=textt(finalindex,1);

graph_overlay_allen_paired(outputpth, spon_pupilhigh_cat_degree,spon_pupillow_cat_degree,'spon_pupilhigh_pupillow','degree_centrality','pupilhigh vs pupillow degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth, spon_pupilhigh_cat_closeness,spon_pupillow_cat_closeness,'spon_pupilhigh_pupillow','closeness_centrality','pupilhigh vs pupillow closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth, spon_pupilhigh_cat_pagerank,spon_pupillow_cat_pagerank,'spon_pupilhigh_pupillow','pagerank_centrality','pupilhigh vs pupillow pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth, spon_pupilhigh_cat_eigenvector,spon_pupillow_cat_eigenvector,'spon_pupilhigh_pupillow','eigenvector_centrality','pupilhigh vs pupillow eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth, spon_pupilhigh_cat_betweenness,spon_pupillow_cat_betweenness,'spon_pupilhigh_pupillow','betweenness_centrality','pupilhigh vs pupillow betweenness Centrality (spon)',parcels_names,length(animals));

graph_overlay_allen_notpaired(outputpth, spon_pupilhigh_cat_degree,spon_pupillow_cat_degree,'spon_pupilhigh_pupillow','degree_centrality','pupilhigh vs pupillow degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_pupilhigh_cat_closeness,spon_pupillow_cat_closeness,'spon_pupilhigh_pupillow','closeness_centrality','pupilhigh vs pupillow closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_pupilhigh_cat_pagerank,spon_pupillow_cat_pagerank,'spon_pupilhigh_pupillow','pagerank_centrality','pupilhigh vs pupillow pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_pupilhigh_cat_eigenvector,spon_pupillow_cat_eigenvector,'spon_pupilhigh_pupillow','eigenvector_centrality','pupilhigh vs pupillow eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_pupilhigh_cat_betweenness,spon_pupillow_cat_betweenness,'spon_pupilhigh_pupillow','betweenness_centrality','pupilhigh vs pupillow betweenness Centrality (spon)',parcels_names,length(animals));
