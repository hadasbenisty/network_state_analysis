function plot_3_bars(M, S, statenames)
CondColors = get_3states_colors;
set(gcf,'renderer','Painters')
hold on
for b = 1:3
    bg=bar(b, M(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
ylim([0.5 0.9]);set(gca,'xtick',1:3);set(gca,'xticklabel',statenames)
h = errorbar(1:3,M, S,'LineStyle','none','LineWidth',0.5);title('Acc All Parcels Per State');
h.Color='k';set(h, 'marker', 'none'); 