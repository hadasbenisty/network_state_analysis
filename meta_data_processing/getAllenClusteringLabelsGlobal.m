function [roiLabelsbyAllen, maskByAllen, regionLabel, isLeftLabel] = getAllenClusteringLabelsGlobal(datapath, bestfile, parcells, maskAllenFrontal, xilinmask)
[R,C] = size(parcells.indicators(:,:,1));
maskByAllen=zeros(R,C);
roiLabelsbyAllen=[];
if isempty(bestfile)
    for a_roi = 1:length(parcells.names)
        maskByAllen(parcells.indicators(:,:,a_roi)==1) = a_roi;
    end
    return;
end
clstrs=unique(maskAllenFrontal);
dat= load(fullfile(datapath, bestfile));


cl = unique(maskAllenFrontal(:));
roilist = setdiff(unique(xilinmask(:)),0);

for roii = 1:size(xilinmask,3)
    P=zeros(R,C);
    pixel_list = xilinmask(:,:,roii) == 1;
    P(pixel_list)=1;
    temp=[];
    for a_roi = 1:length(parcells.names)
        temp(a_roi) = sum(sum(P&parcells.indicators(:,:,a_roi)))/sum(sum(parcells.indicators(:,:,a_roi)));
    end
    [~, roiLabelsbyAllen(roii)] = max(temp);
    maskByAllen(pixel_list) = roiLabelsbyAllen(roii);
    temp=[];
    for a_roi = 1:length(cl)
        temp(a_roi) = sum(sum(P&(maskAllenFrontal==cl(a_roi))))/sum(sum(maskAllenFrontal==cl(a_roi)));
    end
    b=hist(maskAllenFrontal(maskByAllen==roiLabelsbyAllen(roii)),clstrs);
    
    [~, ii] = max(b);
    regionLabel(roii) = clstrs(ii);
    [a,b]=ind2sub([256,256],mean(find(pixel_list)));
    isLeftLabel(roii) = b <= 128;
end
