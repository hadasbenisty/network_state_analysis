function [parcels_names, parcels_region_labels, finalindex, region_lut, allen_map_final_index, parcels_namesall] = get_allen_meta_parcels


isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_region_labels_bilateral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels_bilateral(finalindex);
try
T=readtable('..\meta_data_processing\subregion_list.csv');
parcels_namesall = T.Label;
parcels_names=T.Label(finalindex);
catch
    parcels_names=[];
end
region_lut = {'visual','parietal','temporal','auditory','rs','somato','motor'};
if ~exist('allen_map_final_index.mat','file')
    parcelsallen=load(fullfile('..\meta_data_processing\parcells_updated121519.mat'));
    
    allen_map_final_index = zeros(size(parcelsallen.parcells_new.CombinedParcells));
    for k=1:length(finalindex)
        allen_map_final_index(parcelsallen.parcells_new.indicators(:,:,finalindex(k))==1) = k;
    end
    save('allen_map_final_index', 'allen_map_final_index');
else
    load('allen_map_final_index','allen_map_final_index');
end
