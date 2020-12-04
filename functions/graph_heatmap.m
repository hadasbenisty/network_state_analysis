function graph_heatmap(brain_mask,parcels,condition1,condition2,name,plotitle)
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
%upperlim=max(max(max(condition1)),max(max(condition2)));
difference_cond1_cond2=condition1-condition2;
meandiff=mean(difference_cond1_cond2,2).';
%upperlim=max(meandiff(:)); 
%lowerlim=(min(min(meandiff(:)))); %(min(min(meandiff(:)))-0.00000001);
if contains(name,'eigenvector')&&0
    lowerlim=-0.06;
    upperlim=0.06;
elseif contains(name,'svm')||contains(name,'SVM')
    lowerlim=-0.05;
    upperlim=0.05;
else
    upperlim=max(meandiff(:)); 
    lowerlim=min(meandiff(:));
end
removedindx=setdiff(1:56,finalindex);
maskacc = zeros(size(brain_mask));
maskacc(brain_mask) = 0;
for k = 1:length(finalindex)
    maskacc(parcels(:,:,finalindex(k))==1) = meandiff(:,k);
end
for k = 1:length(removedindx)
    maskacc(parcels(:,:,removedindx(k))==1) = 0;
end
maskacc(:,1:130)=0;

%clims = [lowerlim upperlim];
figure;imagesc(maskacc);
set(gcf,'renderer','painters');
myColorMap = colormap(redblue);
%myColorMap(1,:) = 1;
colormap(myColorMap);
h=colorbar;
caxis([lowerlim upperlim]);
%colormap(fireice);h=colorbar;
ylabel(h, 'Difference in Node Centrality');title(plotitle);axis off
hold on
plot_parcellation_boundaries(parcels(:,:,finalindex));
end

