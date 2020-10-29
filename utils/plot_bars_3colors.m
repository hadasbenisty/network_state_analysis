function plot_bars_3colors(M, S, legstr, xlbels)

figure;
set(gcf,'renderer','Painters')
h=barwitherr(S,M);
CondColors=[0,0,0;0.9961,0.5469,0;0.6953,0.1328,0.1328;0.9961,0.8398,0];
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
legend(legstr)

set(gca,'XTick', 1:length(xlbels));
set(gca,'XTickLabel', xlbels);

set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
set(gca,'XTickLabelRotation',45);
set(gcf,'Position',[1          41        1920         963])