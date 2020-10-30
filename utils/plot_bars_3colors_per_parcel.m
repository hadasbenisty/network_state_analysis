function plot_bars_3colors_per_parcel(M, S, statenames, parcels_names)

figure;
set(gcf,'renderer','Painters')

h=barwitherr(S,M);
CondColors=get_3states_colors;
h(1).EdgeColor = 'none';
h(2).EdgeColor = 'none';
h(3).EdgeColor = 'none';
set(h(1),'FaceColor',CondColors(1,:));
set(h(2),'FaceColor',CondColors(2,:));
set(h(3),'FaceColor',CondColors(3,:));
set(h(1),'FaceAlpha',0.8);
set(h(2),'FaceAlpha',0.8);
set(h(3),'FaceAlpha',0.8);
set(gcf, 'Position',  [150,150, 1500,700]);
legend(statenames)
ylim([0.5 0.9]);
set(gca,'XTick', 1:length(parcels_names));
set(gca,'XTickLabel', parcels_names);

set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
set(gca,'XTickLabelRotation',45);
set(gcf,'Position',[1          41        1500         700])