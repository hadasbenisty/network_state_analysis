function [roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal)

%files = dir(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '*imaging_time_traces_global.mat']);
%load(fullfile(files(1).folder, files(1).name), 'roiLabelsbyAllen', 'regionLabel', 'maskByAllen');
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\gal\' animal '_finalparcellation_Gal.mat'], 'final_par_mask');
% roiLabelsbyAllen_gal = roiLabelsbyAllen.Gal;
% regionLabel_gal = regionLabel.Gal;
% maskByAllen_gal = maskByAllen.Gal;
% maskByGal = zeros(size(maskByAllen.Gal,1));
% for k=1:size(final_par_mask,3)
%     maskByGal(final_par_mask(:,:,k) == 1) = k;
% end


for i = 1:size(final_par_mask,3)
    roiinds = find(final_par_mask(:,:,i) == 1);
            
   
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


