function plot_vals_heatmap(maskacc, clrbarttl , allen_indicators, lowerlim, upperlim, whitenans,colormapvals)
if ~exist('upperlim','var')
upperlim=max(maskacc(:));
lowerlim=min(maskacc(:));
end
if ~exist('domidval','var')
    domidval=true;
end

%clims = [lowerlim upperlim];
imshow(maskacc);
b = imagesc(maskacc);
set(b,'AlphaData',~isnan(maskacc))
 
set(gcf,'renderer','painters');
%myColorMap(1,:) = 1;
colormap(colormapvals);
h=colorbar;
caxis([lowerlim upperlim]);
%colormap(fireice);h=colorbar;
ylabel(h, clrbarttl);axis off
hold on
if ~isempty(allen_indicators)
plot_parcellation_boundaries(allen_indicators);
end
end

