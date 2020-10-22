function [averaged_condition1_eigenvector,averaged_condition2_eigenvector] = permuted_centrality(cond1,cond2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get permuted distribution of centrality measures based on condition
% input:
%   cond1               - name of condition 1 as it is saved in 'X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Permuted_Networks\NetworkAnalysis\'
%   cond2              - name of condition 2
%
% output:
%   averaged_condition_differences_eigenvector          - parcels * animals matrix 
%containing the average hub centrality difference of fake condition 1 and fake condition 2 for each animal and parcel.
%
%   averaged_condition_differences_closeness    - parcels * animals matrix 
%containing the average closeness centrality difference of fake condition 1 and fake condition 2 for each animal and parcel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partial corr
animals={'xw','xx','xz','xs','xt','xu'};
averaged_condition1_eigenvector=[];averaged_condition2_eigenvector=[];
%for each animal, get the average difference across 20 "fake/permtued"
%condition 1 and 2outputfolder=fullfile('X:\Lav\ProcessingDirectory\',animal,'\');
for i=1:length(animals)
    animal=char(animals(i)); 
    k_condition1=[];k_condition2=[];
    for k=1:100
        condition1=strcat(num2str(k),'network_analysis_corr',cond1);condition2=strcat(num2str(k),'network_analysis_corr',cond2);
        %load "fake" condition 1 and condition 2's centrality measures
        condition1_centrality=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\permute\',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
        condition2_centrality=load(strcat('X:\Lav\ProcessingDirectory\',animal,'\permute\',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
        %get fake differences and concatenate them
        k_condition1=cat(2, k_condition1,condition1_centrality.cent_corr.eigenvector); 
        k_condition2=cat(2, k_condition2,condition2_centrality.cent_corr.eigenvector); 
        clearvars condition1_centrality condition2_centrality
    end    
    %for each animal, average over 20 different fake differences to get a
    %mean fake difference.
    averaged_condition1_eigenvector=cat(2, averaged_condition1_eigenvector,nanmean(k_condition1,2)); 
    averaged_condition2_eigenvector=cat(2, averaged_condition2_eigenvector,nanmean(k_condition2,2)); 
    clearvars k_difference_conditions_closeness k_difference_conditions_eigenvector
end
end


