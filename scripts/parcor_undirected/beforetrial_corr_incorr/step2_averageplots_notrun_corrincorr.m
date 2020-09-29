function step2_averageplots_notrun_corrincorr(inputpth, oututpth)
%% corr
clear;close all;
cd('X:\Lav\ProcessingDirectory\parcor_undirected\')
addpath(genpath('X:\Lav\network_state_analysis\functions\'))
animals={'xw','xx','xz','xs','xt','xu'};
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
trial_notrun_correct_cat_degree=[];trial_notrun_incorrect_cat_degree=[];
trial_notrun_correct_cat_closeness=[];trial_notrun_incorrect_cat_closeness=[];
trial_notrun_correct_cat_betweenness=[];trial_notrun_incorrect_cat_betweenness=[];
trial_notrun_correct_cat_pagerank=[];trial_notrun_incorrect_cat_pagerank=[];
trial_notrun_correct_cat_eigenvector=[];trial_notrun_incorrect_cat_eigenvector=[];
condition1='trial_notrun_correct';
condition2='trial_notrun_incorrect';
parcels_region_labels_bitrial_notrun_incorrectral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bitrial_notrun_incorrectral(finalindex);
counter=1;
for i=1:length(animals)
    animal=char(animals(i)); 
    trial_notrun_correct=load(strcat(inputpth,animal,'\','network_analysis_corr',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
    trial_notrun_incorrect=load(strcat(inputpth,animal,'\','network_analysis_corr',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
    trial_notrun_correct_cat_degree=cat(2, trial_notrun_correct_cat_degree,trial_notrun_correct.cent_corr.degree); 
    trial_notrun_correct_cat_closeness=cat(2, trial_notrun_correct_cat_closeness,trial_notrun_correct.cent_corr.closeness); 
    trial_notrun_correct_cat_betweenness=cat(2, trial_notrun_correct_cat_betweenness,trial_notrun_correct.cent_corr.betweenness);    
    trial_notrun_correct_cat_pagerank=cat(2, trial_notrun_correct_cat_pagerank,trial_notrun_correct.cent_corr.pagerank); 
    trial_notrun_correct_cat_eigenvector=cat(2, trial_notrun_correct_cat_eigenvector,trial_notrun_correct.cent_corr.eigenvector); 
     subplot(2,length(animals),counter)
    p=plot(trial_notrun_correct.G_corr);title(strcat(condition1,animal));counter=counter+1;
    p.NodeCData=parcels_region_labels;
    p.NodeLabel = {};
    subplot(2,length(animals),counter)
    n=plot(trial_notrun_incorrect.G_corr);title(strcat(condition2,animal));
    n.NodeCData=parcels_region_labels;counter=counter+1;
    n.NodeLabel = {};
    %mysave(gcf,strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_network\',animal,'\','trial_notrun_correct_trial_notrun_incorrect_graph'),'all')
    
    trial_notrun_incorrect_cat_degree=cat(2, trial_notrun_incorrect_cat_degree,trial_notrun_incorrect.cent_corr.degree); 
    trial_notrun_incorrect_cat_closeness=cat(2, trial_notrun_incorrect_cat_closeness,trial_notrun_incorrect.cent_corr.closeness); 
    trial_notrun_incorrect_cat_betweenness=cat(2, trial_notrun_incorrect_cat_betweenness,trial_notrun_incorrect.cent_corr.betweenness);       
    trial_notrun_incorrect_cat_pagerank=cat(2, trial_notrun_incorrect_cat_pagerank,trial_notrun_incorrect.cent_corr.pagerank); 
    trial_notrun_incorrect_cat_eigenvector=cat(2, trial_notrun_incorrect_cat_eigenvector,trial_notrun_incorrect.cent_corr.eigenvector); 
    clearvars trial_notrun_correct trial_notrun_incorrect
end
[~,textt]=xlsread('X:\Hadas\Meso-imaging\Antara\preprocessing\meso_processing-master\parcellation\AllenParcellationLan\allanParcellationTiffs\subregion_list.csv');
textt(1,:)=[];
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_names=textt(finalindex,1);
graph_overlay_allen_paired(oututpth, trial_notrun_correct_cat_pagerank,trial_notrun_incorrect_cat_pagerank,'trial_notrun_correct_trial_notrun_incorrect','pagerank_centrality','trial_notrun_correct vs trial_notrun_incorrect pagerank Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_notpaired(oututpth, trial_notrun_correct_cat_eigenvector,trial_notrun_incorrect_cat_eigenvector,'trial_notrun_correct_trial_notrun_incorrect','eigenvector_centrality','trial_notrun_correct vs trial_notrun_incorrect eigenvector Centrality (spon)',parcels_names,length(animals));

graph_overlay_allen_paired(oututpth, trial_notrun_correct_cat_degree,trial_notrun_incorrect_cat_degree,'trial_notrun_correct_trial_notrun_incorrect','degree_centrality','trial_notrun_correct vs trial_notrun_incorrect degree Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(oututpth, trial_notrun_correct_cat_closeness,trial_notrun_incorrect_cat_closeness,'trial_notrun_correct_trial_notrun_incorrect','closeness_centrality','trial_notrun_correct vs trial_notrun_incorrect closeness Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(oututpth, trial_notrun_correct_cat_eigenvector,trial_notrun_incorrect_cat_eigenvector,'trial_notrun_correct_trial_notrun_incorrect','eigenvector_centrality','trial_notrun_correct vs trial_notrun_incorrect eigenvector Centrality (spon)',parcels_names,length(animals));
graph_overlay_allen_paired(oututpth, trial_notrun_correct_cat_betweenness,trial_notrun_incorrect_cat_betweenness,'trial_notrun_correct_trial_notrun_incorrect','betweenness_centrality','trial_notrun_correct vs trial_notrun_incorrect betweenness Centrality (spon)',parcels_names,length(animals));

