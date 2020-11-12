function [perm1_cat,perm2_cat] = load_permutedplots(cond1,cond2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get permuted distribution of centrality measures based on condition
% input:
%   cond1       - name of condition 1 
%   cond2       - name of condition 2
%
% output:
%   perm1_cat   - parcels * animals matrix 
%containing the average eigenvector centrality of fake condition 1
%
%   perm2_cat    - parcels * animals matrix 
%containing the average eigenvector centrality of fake condition 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
animals={'xw','xx','xz','xs','xt','xu'};
perm1_cat=[];perm2_cat=[];
% spon_run_cat_authorities=[];spon_notrun_cat_authorities=[];
for i=1:length(animals)
    animal=char(animals(i));
    k_cond1=[];k_cond2=[];
    for k=1:100
        condition1=strcat(num2str(k),'network_analysis_corr_perm',cond1);condition2=strcat(num2str(k),'network_analysis_corr_perm',cond2);
        
        cond1data=load(strcat('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
        cond2data=load(strcat('X:\Lav\ProcessingDirectory_Oct2020\',animal,'\',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');  
        k_cond1=cat(2, k_cond1,cond1data.cent_corr.eigenvector);
        k_cond2=cat(2, k_cond2,cond2data.cent_corr.eigenvector);
    end    
    perm1_cat=cat(2, perm1_cat,nanmean(k_cond1,2));
    perm2_cat=cat(2, perm2_cat,nanmean(k_spon_notrun_cat_hubs,2));
    clearvars cond1data cond2data
end
end