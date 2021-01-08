function graph_overlay_allen_2conditions(strsuffix, outputpth,cond1_cat,cond3_cat,condition,name,plottitle,parcels_names,legstr)
    %calculate average and error of condition 1
    num = sum(~isnan(cond1_cat(1,:)));
    condition1_meanvec_parcel=nanmean(cond1_cat,2);
    condition1_sevec_parcel=nanstd(cond1_cat,0,2)./sqrt(num-1);
    
    
    condition3_meanvec_parcel=nanmean(cond3_cat,2);
    condition3_sevec_parcel=nanstd(cond3_cat,0,2)./sqrt(num-1);
    CondColors=[0,0,0;0.9961,0.5469,0;1,0,0];
    CondColors=CondColors([1 3],:);
    %get colors
    figure;
    set(gcf,'renderer','Painters')
    y = cat(2,condition1_meanvec_parcel,condition3_meanvec_parcel);         % random y values (3 groups of 4 parameters)
    errY = cat(2,condition1_sevec_parcel,condition3_sevec_parcel);          % 10% error
    h = barwitherr(errY, y);% Plot with errorbars
    h(1).EdgeColor = 'none';
    h(2).EdgeColor = 'none';
    set(h(1),'FaceColor',CondColors(1,:));
    set(h(2),'FaceColor',CondColors(2,:));
    set(h(1),'FaceAlpha',0.8);
    set(h(2),'FaceAlpha',0.8);
    set(gcf, 'Position',  [150,150, 1500,700]);
   
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    ylabel('Node Centrality')

    set(gcf, 'Position',  [150,150, 2000,500]);
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    ylabel('Node Centrality');title(plottitle);
    legend(legstr);
    mkNewDir(fullfile(outputpth,condition));
    mysave(gcf, fullfile(outputpth,condition,strcat(name,'2conditions', strsuffix)), 'all');
end



