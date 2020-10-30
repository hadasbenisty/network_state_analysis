function [roiLabelsbyAllen_gal, regionLabel_gal, maskByAllen_gal, maskByGal] = get_gal_parcels_lables(animal)

files = dir(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '*imaging_time_traces_global.mat']);
load(fullfile(files(1).folder, files(1).name), 'roiLabelsbyAllen', 'regionLabel', 'maskByAllen');
load(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\gal\' animal '_finalparcellation_Gal.mat'], 'final_par_mask');
roiLabelsbyAllen_gal = roiLabelsbyAllen.Gal;
regionLabel_gal = regionLabel.Gal;
maskByAllen_gal = maskByAllen.Gal;
maskByGal = zeros(size(maskByAllen.Gal,1));
for k=1:size(final_par_mask,3)
    maskByGal(final_par_mask(:,:,k) == 1) = k;
end