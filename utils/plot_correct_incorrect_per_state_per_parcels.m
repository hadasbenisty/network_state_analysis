function ha = plot_correct_incorrect_per_state_per_parcels(M, S, parcels_names, statenames)
CondColors=get_3states_colors;
figure;
for state_i = 1:length(statenames)
    ha(state_i) = subplot(3,1,state_i);
    
    set(gcf,'renderer','Painters');
    h=barwitherr(S(:,:, state_i),M(:,:, state_i));title(statenames{state_i});
    h(1).EdgeColor = 'none';
    h(2).EdgeColor = 'none';
    set(h(1),'FaceColor',CondColors(state_i,:));
    set(h(2),'FaceColor',CondColors(state_i,:)); 
    set(h(1),'FaceAlpha',1);
    set(h(2),'FaceAlpha',0.4); 
    set(gca,'xtick',1:23)
    %set(gcf, 'Position',  [1,1, 700,1000]); 
    set(gcf, 'Position',  [1,41, 1920,963]);
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',10)
    set(gca,'XTickLabelRotation',45);
   
end