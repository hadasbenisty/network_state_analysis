function step2_averageplots(inputpth, outputpth) 
%% corr
% clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
spon_facehigh_cat_degree=[];spon_facelow_cat_degree=[];
spon_facehigh_cat_closeness=[];spon_facelow_cat_closeness=[];
spon_facehigh_cat_betweenness=[];spon_facelow_cat_betweenness=[];
spon_facehigh_cat_pagerank=[];spon_facelow_cat_pagerank=[];
spon_facehigh_cat_eigenvector=[];spon_facelow_cat_eigenvector=[];
condition1='facehigh';condition2='facelow';
parcels_region_labels_bifacelowral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bifacelowral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    spon_facehigh=load(strcat(inputpth, animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    spon_facelow=load(strcat(inputpth,animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    spon_facehigh_cat_degree=cat(2, spon_facehigh_cat_degree,spon_facehigh.cent_corr.degree); 
    spon_facehigh_cat_closeness=cat(2, spon_facehigh_cat_closeness,spon_facehigh.cent_corr.closeness); 
    spon_facehigh_cat_betweenness=cat(2, spon_facehigh_cat_betweenness,spon_facehigh.cent_corr.betweenness);    
    spon_facehigh_cat_pagerank=cat(2, spon_facehigh_cat_pagerank,spon_facehigh.cent_corr.pagerank); 
    spon_facehigh_cat_eigenvector=cat(2, spon_facehigh_cat_eigenvector,spon_facehigh.cent_corr.eigenvector); 
    subplot(2,length(animals),counter)
    p=plot(spon_facehigh.G_corr);title(strcat(condition1,animal));counter=counter+1;
    p.NodeCData=parcels_region_labels;
    p.NodeLabel = {};
    subplot(2,length(animals),counter)
    n=plot(spon_facelow.G_corr);title(strcat(condition2,animal));
    n.NodeCData=parcels_region_labels;counter=counter+1;
    n.NodeLabel = {};
    %mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_spon_network\',animal,'\','facehigh_facelow_graph'),'all')
    
    spon_facelow_cat_degree=cat(2, spon_facelow_cat_degree,spon_facelow.cent_corr.degree); 
    spon_facelow_cat_closeness=cat(2, spon_facelow_cat_closeness,spon_facelow.cent_corr.closeness); 
    spon_facelow_cat_betweenness=cat(2, spon_facelow_cat_betweenness,spon_facelow.cent_corr.betweenness);       
    spon_facelow_cat_pagerank=cat(2, spon_facelow_cat_pagerank,spon_facelow.cent_corr.pagerank); 
    spon_facelow_cat_eigenvector=cat(2, spon_facelow_cat_eigenvector,spon_facelow.cent_corr.eigenvector); 
    diffusionmap_twoconditions(spon_facehigh.W_corr,spon_facelow.W_corr,spon_facehigh.cent_corr.eigenvector,spon_facelow.cent_corr.eigenvector,condition1,condition2,animal,'eigenvector')

    clearvars spon_facehigh spon_facelow
end
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_names=textt(finalindex,1);
mkNewDir(outputpth);
%[avg_cond1,avg_cond2] = permuted_centrality('highpup','lowpup');
%graph_overlay_allen_paired_permuted(avg_cond1,avg_cond2,outputpth, spon_facehigh_cat_eigenvector,spon_facelow_cat_eigenvector,'spon_facehigh_facelow','eigenvector_centrality','facehigh vs facelow eigenvector Centrality (spon)',parcels_names,length(animals));

% graph_overlay_allen_paired(outputpth, spon_facehigh_cat_degree,spon_facelow_cat_degree,'spon_facehigh_facelow','degree_centrality','facehigh vs facelow degree Centrality (spon)',parcels_names,length(animals));
% graph_overlay_allen_paired(outputpth, spon_facehigh_cat_closeness,spon_facelow_cat_closeness,'spon_facehigh_facelow','closeness_centrality','facehigh vs facelow closeness Centrality (spon)',parcels_names,length(animals));
% graph_overlay_allen_paired(outputpth, spon_facehigh_cat_pagerank,spon_facelow_cat_pagerank,'spon_facehigh_facelow','pagerank_centrality','facehigh vs facelow pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(outputpth, spon_facehigh_cat_eigenvector,spon_facelow_cat_eigenvector,'spon_facehigh_facelow','eigenvector_centrality','facehigh vs facelow eigenvector Centrality (spon)',parcels_names,length(animals));
% graph_overlay_allen_paired(outputpth, spon_facehigh_cat_betweenness,spon_facelow_cat_betweenness,'spon_facehigh_facelow','betweenness_centrality','facehigh vs facelow betweenness Centrality (spon)',parcels_names,length(animals));

braininfo=load('X:\Lav\network_state_analysis\utils\brain_mask.mat');
parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');

graph_heatmap(outputpth,braininfo.brain_mask,parcelsallen.parcells_new.indicators, spon_facehigh_cat_degree,spon_facelow_cat_degree,'spon_facehigh_facelow','degree_centrality','facehigh vs facelow degree Centrality (spon)');
graph_heatmap(outputpth,braininfo.brain_mask,parcelsallen.parcells_new.indicators, spon_facehigh_cat_closeness,spon_facelow_cat_closeness,'spon_facehigh_facelow','closeness_centrality','facehigh vs facelow closeness Centrality (spon)');
graph_heatmap(outputpth,braininfo.brain_mask,parcelsallen.parcells_new.indicators, spon_facehigh_cat_pagerank,spon_facelow_cat_pagerank,'spon_facehigh_facelow','pagerank_centrality','facehigh vs facelow pagerank Centrality (spon)');
graph_heatmap(outputpth,braininfo.brain_mask,parcelsallen.parcells_new.indicators, spon_facehigh_cat_eigenvector,spon_facelow_cat_eigenvector,'spon_facehigh_facelow','eigenvector_centrality','facehigh vs facelow eigenvector Centrality (spon)');
graph_heatmap(outputpth,braininfo.brain_mask,parcelsallen.parcells_new.indicators, spon_facehigh_cat_betweenness,spon_facelow_cat_betweenness,'spon_facehigh_facelow','betweenness_centrality','facehigh vs facelow betweenness Centrality (spon)');

graph_overlay_allen_notpaired(outputpth, spon_facehigh_cat_degree,spon_facelow_cat_degree,'spon_facehigh_facelow','degree_centrality','facehigh vs facelow degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_facehigh_cat_closeness,spon_facelow_cat_closeness,'spon_facehigh_facelow','closeness_centrality','facehigh vs facelow closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_facehigh_cat_pagerank,spon_facelow_cat_pagerank,'spon_facehigh_facelow','pagerank_centrality','facehigh vs facelow pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_facehigh_cat_eigenvector,spon_facelow_cat_eigenvector,'spon_facehigh_facelow','eigenvector_centrality','facehigh vs facelow eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(outputpth, spon_facehigh_cat_betweenness,spon_facelow_cat_betweenness,'spon_facehigh_facelow','betweenness_centrality','facehigh vs facelow betweenness Centrality (spon)',parcels_names,length(animals));

