function P = scores_to_heatmap_allen(Xvec)



[parcels_names, ~, ~, ~, allen_map_final_index] = get_allen_meta_parcels;


P=nan(256);
for parcel_i = 1:length(parcels_names)
    P(allen_map_final_index==parcel_i) = Xvec(parcel_i);
end


      