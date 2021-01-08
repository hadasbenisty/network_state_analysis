function plot_2_bars(M, S, statenames)
CondColors = get_3states_colors;
CondColors=CondColors([1 3], :);
set(gcf,'renderer','Painters')
hold on
for b = 1:2
    bg=bar(b, M(b), 'FaceColor',  CondColors(b,:), 'EdgeColor', 'none', 'BarWidth', 0.6);hold on;
    bg.FaceAlpha = 0.8;
end
set(gca,'xtick',1:2);set(gca,'xticklabel',statenames)
h = errorbar(1:2,M, S,'LineStyle','none','LineWidth',0.5);
h.Color='k';set(h, 'marker', 'none'); 