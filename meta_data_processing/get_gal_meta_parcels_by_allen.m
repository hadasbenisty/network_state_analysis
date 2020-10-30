function [parcels_names_gal, finalindex_gal, maskByGalnew, regionLabel_galnew] = get_gal_meta_parcels_by_allen(parcels_names, finalindex,...
    roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal)


finalindex_gal = [];
parcels_names_gal = [];
regionLabel_galnew=[];
maskByGalnew = zeros(size(maskByAllen_gal));
l=1;
for i=1:length(finalindex)
    inds = find(roiLabelsbyAllen_gal == finalindex(i));
   finalindex_gal = cat(2, finalindex_gal, inds);  
   regionLabel_galnew = cat(2, regionLabel_galnew, regionLabel_gal(inds));
   for k=1:length(inds)
       maskByGalnew(maskByGal==inds(k))=l;
       l=l+1;
       parcels_names_gal{end+1} = parcels_names{i};
   end
end

