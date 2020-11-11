function plot_heatmap(A, lowlim, uplim, ttl, parcels)

A(isnan(A)) = 0;
    imagesc(A,[-0.04 0.04]);
set(gcf,'renderer','painters');
myColorMap = colormap(redblue);
%myColorMap(1,:) = 1;
colormap(myColorMap);
h=colorbar;ylabel(h, 'Difference in Node Centrality');
caxis([lowlim uplim]);
title(ttl);axis off
hold on;axis off;
plot_parcellation_boundaries(parcels);