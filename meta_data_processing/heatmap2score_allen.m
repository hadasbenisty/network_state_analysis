function scores = heatmap2score_allen(heatmap)


[parcels_names, ~, ~, ~, allen_map_final_index] = get_allen_meta_parcels;
scores=nan(length(parcels_names),size(heatmap,3));

for parcel_i = 1:length(parcels_names)
    inds = allen_map_final_index==parcel_i;
    for k=1:size(heatmap,3)
        P=heatmap(:,:,k);
        scores(parcel_i,k) = nanmean(P(inds));
    end
end

