function step2_averageplots_3states(inputpth, outputpth) 
%% corr
% clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
addpath(genpath('X:\Lav\network_state_analysis\utils'))

animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
spon_pupilhigh_cat_degree=[];spon_pupillow_cat_degree=[];
spon_pupilhigh_cat_closeness=[];spon_pupillow_cat_closeness=[];
spon_pupilhigh_cat_betweenness=[];spon_pupillow_cat_betweenness=[];
spon_pupilhigh_cat_pagerank=[];spon_pupillow_cat_pagerank=[];
spon_pupilhigh_cat_eigenvector=[];spon_pupillow_cat_eigenvector=[];
spon_run_cat_degree=[];
spon_run_cat_closeness=[];
spon_run_cat_betweenness=[];
spon_run_cat_pagerank=[];
spon_run_cat_eigenvector=[];
condition1='pupilhigh';condition2='pupillow';condition3='run';
parcels_region_labels_bipupillowral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bipupillowral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    spon_pupilhigh=load(strcat(inputpth, animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    spon_pupillow=load(strcat(inputpth,animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    spon_run=load(strcat(inputpth,animal,'\','network_analysis_corr',condition3),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
%% state 1
    spon_run_cat_degree=cat(2, spon_run_cat_degree,spon_run.cent_corr.degree); 
    spon_run_cat_closeness=cat(2, spon_run_cat_closeness,spon_run.cent_corr.closeness); 
    spon_run_cat_betweenness=cat(2, spon_run_cat_betweenness,spon_run.cent_corr.betweenness);    
    spon_run_cat_pagerank=cat(2, spon_run_cat_pagerank,spon_run.cent_corr.pagerank); 
    spon_run_cat_eigenvector=cat(2, spon_run_cat_eigenvector,spon_run.cent_corr.eigenvector); 
%% state 2
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
%% state 3    
    spon_pupillow_cat_degree=cat(2, spon_pupillow_cat_degree,spon_pupillow.cent_corr.degree); 
    spon_pupillow_cat_closeness=cat(2, spon_pupillow_cat_closeness,spon_pupillow.cent_corr.closeness); 
    spon_pupillow_cat_betweenness=cat(2, spon_pupillow_cat_betweenness,spon_pupillow.cent_corr.betweenness);       
    spon_pupillow_cat_pagerank=cat(2, spon_pupillow_cat_pagerank,spon_pupillow.cent_corr.pagerank); 
    spon_pupillow_cat_eigenvector=cat(2, spon_pupillow_cat_eigenvector,spon_pupillow.cent_corr.eigenvector); 
    %diffusionmap_twoconditions(spon_pupilhigh.W_corr,spon_pupillow.W_corr,spon_pupilhigh.cent_corr.eigenvector,spon_pupillow.cent_corr.eigenvector,condition1,condition2,animal,'eigenvector')

    clearvars spon_pupilhigh spon_pupillow
end
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_names=textt(finalindex,1);
mkNewDir(outputpth);
%[avg_cond1,avg_cond2] = permuted_centrality('highpup','lowpup');
%graph_overlay_allen_paired_permuted(avg_cond1,avg_cond2,outputpth, spon_pupilhigh_cat_eigenvector,spon_pupillow_cat_eigenvector,'spon_pupilhigh_pupillow','eigenvector_centrality','pupilhigh vs pupillow eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_3conditions(outputpth, spon_pupilhigh_cat_eigenvector,spon_pupillow_cat_eigenvector,spon_run_cat_eigenvector,'spon_threestates','eigenvector_centrality','3 states eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_3conditions(outputpth, spon_pupilhigh_cat_betweenness,spon_pupillow_cat_betweenness,spon_run_cat_betweenness,'spon_threestates','betweenness_centrality','3 states betweenness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_3conditions(outputpth, spon_pupilhigh_cat_degree,spon_pupillow_cat_degree,spon_run_cat_degree,'spon_threestates','degree_centrality','3 states degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_3conditions(outputpth, spon_pupilhigh_cat_closeness,spon_pupillow_cat_closeness,spon_run_cat_closeness,'spon_threestates','closeness_centrality','3 states closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_3conditions(outputpth, spon_pupilhigh_cat_pagerank,spon_pupillow_cat_pagerank,spon_run_cat_pagerank,'spon_threestates','pagerank_centrality','3 states pagerank Centrality (spon)',parcels_names,length(animals));



