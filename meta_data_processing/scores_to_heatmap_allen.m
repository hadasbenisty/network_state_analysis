function Pvec = scores_to_heatmap_allen(Xmat, tonorm)


[parcels_names, ~, ~, ~, allen_map_final_index] = get_allen_meta_parcels;
Pvec=nan(256,256,size(Xmat,2));
for kk=1:size(Xmat,2)
    P=nan(256);
    for parcel_i = 1:length(parcels_names)
        if tonorm
            P(allen_map_final_index==parcel_i) = Xmat(parcel_i,kk)/sqrt(size(Xmat,1));
        else
            P(allen_map_final_index==parcel_i) = Xmat(parcel_i,kk);
        end
    end
    Pvec(:,:,kk)=P;
end


      