function P = scores_to_heatmap_gal(Xvec, animal)



[parcels_names, ~, finalindex] = get_allen_meta_parcels;
[roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal);

[parcels_names_gal, ~, maskByGal] = get_gal_meta_parcels_by_allen(parcels_names, finalindex,...
    roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal);


P=nan(256);
for parcel_i = 1:length(parcels_names_gal)
    P(maskByGal==parcel_i) = Xvec(parcel_i);
end


      