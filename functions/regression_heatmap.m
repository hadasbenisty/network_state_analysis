function regression_heatmap(brain_mask,parcels,betas,plotitle)
%isleftlabel=2:2:56;
%finalindex=setdiff(1:56,[21:26 53:56]);
%finalindex=intersect(isleftlabel,toremove);
meandiff=nanmean(betas,2).';
removedindx=[21:26 53:56];
maskacc = zeros(size(brain_mask));
maskacc(brain_mask) = 0;

for k = 1:56
    maskacc(parcels(:,:,k)==1) = meandiff(k);
end
for k = 1:length(removedindx)
    maskacc(parcels(:,:,removedindx(k))==1) = 0;
end
%maskacc(:,1:130)=0;

figure;imagesc(maskacc);
set(gcf,'renderer','painters');
myColorMap = colormap(redblue);
colormap(myColorMap);
h=colorbar;
%caxis([lowerlim upperlim]);
ylabel(h, 'R^2');title(plotitle);axis off
hold on
plot_parcellation_boundaries(parcels(:,:,:));
end

