function [parcels_names_grid, parcels_region_labels_grid, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, G)

[parcels_names, parcels_region_labels, finalindex, region_lut, allen_map_final_index] = get_allen_meta_parcels;
grid_map_final_index=zeros(size(allen_map_final_index));
final_index_grid = [];parcels_names_grid = [];parcels_region_labels_grid=[];labelsbyallen=[];
for i = 1:size(par_inds,2)
    [II,JJ]=meshgrid(par_inds(1,i):par_inds(1,i)+G-1,par_inds(2,i):par_inds(2,i)+G-1);
            roiinds = sub2ind([256 256],II(:),JJ(:));
            
   
    if any(allen_map_final_index(roiinds)~=0)
        votes = allen_map_final_index(roiinds);
        votes=votes(votes~=0);
        [votes, bins] = hist(votes, unique(votes));
        [~,maxind] = max(votes);
        selvote = bins(maxind);
        final_index_grid(end+1) = i;
        if ~isempty(parcels_names)
        parcels_names_grid{end+1} = parcels_names{selvote};
        end
        grid_map_final_index(roiinds) = i;
        labelsbyallen(end+1) = selvote;
        parcels_region_labels_grid(end+1) = parcels_region_labels(selvote);
    end
end
