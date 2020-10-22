function [averaged_condition_differences_hubs,averaged_condition_differences_authorities] = load_permutedplots(cond1,cond2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get permuted distribution of centrality measures based on condition
% input:
%   cond1               - name of condition 1 as it is saved in 'X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Permuted_Networks\NetworkAnalysis\'
%   cond2              - name of condition 2
%
% output:
%   averaged_condition_differences_hubs          - parcels * animals matrix 
%containing the average hub centrality difference of fake condition 1 and fake condition 2 for each animal and parcel.
%
%   averaged_condition_differences_authorities    - parcels * animals matrix 
%containing the average authorities centrality difference of fake condition 1 and fake condition 2 for each animal and parcel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partial corr
animals={'xw','xx','xz','xs','xt','xu'};
averaged_condition_differences_hubs=[];averaged_condition_differences_authorities=[];
%for each animal, get the average difference across 20 "fake/permtued"
%condition 1 and 2
for i=1:length(animals)
    animal=char(animals(i)); 
    k_difference_conditions_authorities=[];k_difference_conditions_hubs=[];
    for k=1:20
        condition1=strcat(num2str(k),'network_analysis_corr',cond1);condition2=strcat(num2str(k),'network_analysis_corr',cond2);
        %load "fake" condition 1 and condition 2's centrality measures
        condition1_centrality=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Permuted_Networks\NetworkAnalysis\',animal,'\',condition1),'W_corr','indic_corr','cent_corr','G_corr','names_corr');
        condition2_centrality=load(strcat('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Permuted_Networks\NetworkAnalysis\',animal,'\',condition2),'W_corr','indic_corr','cent_corr','G_corr','names_corr');   
        %get fake differences and concatenate them
        k_difference_conditions_authorities=cat(2, k_difference_conditions_authorities,(condition1_centrality.cent_corr.authorities-condition2_centrality.cent_corr.authorities)); 
        k_difference_conditions_hubs=cat(2, k_difference_conditions_hubs,(condition1_centrality.cent_corr.hubs-condition2_centrality.cent_corr.hubs)); 
        clearvars condition1_centrality condition2_centrality
    end    
    %for each animal, average over 20 different fake differences to get a
    %mean fake difference.
    averaged_condition_differences_hubs=cat(2, averaged_condition_differences_hubs,nanmean(k_difference_conditions_hubs,2)); 
    averaged_condition_differences_authorities=cat(2, averaged_condition_differences_authorities,nanmean(k_difference_conditions_authorities,2)); 
    clearvars k_difference_conditions_authorities k_difference_conditions_hubs
end
end

