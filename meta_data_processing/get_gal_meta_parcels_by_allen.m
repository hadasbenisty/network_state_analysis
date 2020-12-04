function [parcels_names_gal, finalindex_gal, maskByGal, regionLabel_gal, roiLabelsbyAllen_gal] = get_gal_meta_parcels_by_allen(parcels_names, maskByAllen, ...
    regionLabel, animal)%;roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal)

load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\gal\' animal '_finalparcellation_Gal.mat'], 'final_par_mask');

regionLabel_gal=[];
maskByGal = zeros(size(maskByAllen));
finalindex_gal=[];l=1;parcels_names_gal={};
for i=1:size(final_par_mask, 3)
    inds = find(final_par_mask(:,:,i) >0);
    allenvotes = maskByAllen(inds);
    if all(maskByAllen(inds) == 0)
        continue;
    end
    [b,a]=hist(allenvotes, unique(allenvotes));    
    
   finalindex_gal(l) = i;
   maskByGal(inds) = l;
   [~,mmax] = max(b);
   selparcel_allen = a(mmax);
   parcels_names_gal(l) = parcels_names(selparcel_allen);
   roiLabelsbyAllen_gal(l) = selparcel_allen;
   regionLabel_gal(l) = regionLabel(selparcel_allen);
   l=l+1;
end

