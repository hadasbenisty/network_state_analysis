function plot_bars_3colors_per_parcel(M, S, statenames, parcels_names)

figure;
set(gcf,'renderer','Painters')

h=barwitherr(S,M);
CondColors=get_3states_colors;
for i=1:length(h)
h(i).EdgeColor = 'none';
set(h(i),'FaceColor',CondColors(i,:));
set(h(i),'FaceAlpha',0.8);


end
set(gcf, 'Position',  [150,150, 1500,700]);
legend(statenames)
ylim([0.5 0.9]);
set(gca,'XTick', 1:length(parcels_names));
set(gca,'XTickLabel', parcels_names);

set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
set(gca,'XTickLabelRotation',45);
set(gcf,'Position',[1          41        1500         700])