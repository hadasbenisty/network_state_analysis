function Pmat = scores_to_heatmap_gal(Xvec, animal, tonorm)
if ~exist('tonorm','var')
    tonorm = true;
end
[~, allen_parcels] = getParcellsByLansAllansAtlas;

[parcels_names.Allen, ~, finalindex.Allen, regionLabel.nameslegend, maskByAllen.Allen] = get_allen_meta_parcels;
regionLabel.Allen = allen_parcels.regionNum;
regionLabel.Allen=regionLabel.Allen(finalindex.Allen);
[parcels_names.Gal, ~, maskByAllen.Gal] = get_gal_meta_parcels_by_allen(parcels_names.Allen, maskByAllen.Allen, ...
    regionLabel.Allen, animal);
        
Pmat=nan(256,256, size(Xvec,2));
for kk=1:size(Xvec,2)
P=nan(256);
for parcel_i = 1:length(parcels_names.Gal)
    if tonorm
    P(maskByAllen.Gal==parcel_i) = Xvec(parcel_i)/sqrt(length(Xvec));
    else
        P(maskByAllen.Gal==parcel_i) = Xvec(parcel_i);
    end
end
Pmat(:,:,kk) = P;
end

      