function plot_vals_heatmap(maskacc, clrbarttl , allen_indicators, lowerlim, upperlim, whitenans,colormapvals)
if ~exist('upperlim','var')
upperlim=max(maskacc(:));
lowerlim=min(maskacc(:));
end
if ~exist('domidval','var')
    domidval=true;
end
if ~exist('colormapvals','var')
    colormapvals = colormap(redblue);
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
mid = mean([lowerlim upperlim]);
set(h,'Ticks', [lowerlim mid upperlim])
if exist('clrbarttl','var')
ylabel(h, clrbarttl);
end
axis off
hold on
if exist('allen_indicators','var') &&~isempty(allen_indicators)
plot_parcellation_boundaries(allen_indicators);
end
end

