function graph_overlay_allen_paired_permuted(condition_perm1,condition_perm2,outputpth,condition1,condition2,condition,name,plottitle,parcels_names,num)
isleftlabel=2:2:56;
toremove=setdiff(1:56,[21:26 53:56]);
finalindex=intersect(isleftlabel,toremove);
parcels_region_labels=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
parcels_region_labels=parcels_region_labels(finalindex);
RegionColors=[0 0.4470 0.7410;1 1 0.1294;0.4706 0.6706 0.3020;0.8 1 0;0 0.6 0.9;0.9294 0.6902 0.1294;0.9294 0.4000 0.1294;0.9294 0.1020 0.1294];
%compute real centrality differences
difference_cond1_cond2=condition1-condition2;
%average across animals
meandiff=mean(difference_cond1_cond2,2);
%sediff=std(difference_cond1_cond2,0,2)./sqrt(num-1);

%average permuted centrality differences across animals
difference_perm_cond1_cond2=condition_perm1-condition_perm2;
perm_meandiff=abs(mean(difference_perm_cond1_cond2,2));
%perm_sediff=std(difference_perm_cond1_cond2,0,2)./sqrt(num-1);

figure;
set(gcf,'renderer','Painters')
subplot(1,1,1)
hold on
for b = 1:23
    indx=parcels_region_labels(b);
    bar(b, meandiff(b), 'FaceColor',  RegionColors(indx,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
end
hold on
for b = 1:23
    bar(b, perm_meandiff(b), 'FaceColor',  ([17 17 17]/255), 'FaceAlpha',0.5,'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
end
hold on
for b = 1:23
    bar(b, -(perm_meandiff(b)), 'FaceColor',  ([17 17 17]/255), 'FaceAlpha',0.5,'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
end
view([-90 90]);
% show standard deviation on top
% h = errorbar(1:23,meandiff, sediff,'LineStyle','none','LineWidth',0.5);
ylabel('Difference inNode Centrality');title(plottitle);
% h.Color='k';h.CapSize = 0;
% set(h, 'marker', 'none'); % remove marker
% hold all
% plot(find(htest<=0.05), meandiff(htest<=0.05),'k*', 'MarkerSize',12);hold on;
% h2 = errorbar(1:23,perm_meandiff, perm_sediff,'LineStyle','none','LineWidth',0.5);
% hold on;
view([-90 90]);
if contains(name,'eigenvector')
ylim([-0.06 0.06]);
end
set(gcf, 'Position',  [150,150, 900,600]);
set(gca,'xtick',1:23)
set(gca,'xticklabel',parcels_names)
set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
%set(gca,'XTickLabelRotation',45);

mkdir(condition);
mysave(gcf, fullfile(outputpth,condition,name), 'all');
end


