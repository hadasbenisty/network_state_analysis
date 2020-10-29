function graph_overlay_allen_3conditions(strsuffix, outputpth,cond1_cat,cond2_cat,cond3_cat,condition,name,plottitle,parcels_names,num, legstr)
    %calculate average and error of condition 1
    condition1_meanvec_parcel=mean(cond1_cat,2);
    condition1_sevec_parcel=std(cond1_cat,0,2)./sqrt(num-1);
    %calculate average and error of condition 2
    condition2_meanvec_parcel=mean(cond2_cat,2);
    condition2_sevec_parcel=std(cond2_cat,0,2)./sqrt(num-1);
    
    condition3_meanvec_parcel=mean(cond3_cat,2);
    condition3_sevec_parcel=std(cond3_cat,0,2)./sqrt(num-1);
    CondColors=[0,0,0;0.9961,0.5469,0;0.6953,0.1328,0.1328;0.9961,0.8398,0];
    
    %get colors
    figure;
    set(gcf,'renderer','Painters')
    y = cat(2,condition1_meanvec_parcel,condition2_meanvec_parcel,condition3_meanvec_parcel);         % random y values (3 groups of 4 parameters)
    errY = cat(2,condition1_sevec_parcel,condition2_sevec_parcel,condition3_sevec_parcel);          % 10% error
    h = barwitherr(errY, y);% Plot with errorbars
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
   
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    ylabel('Node Centrality')

    set(gcf, 'Position',  [150,150, 1500,700]);
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    ylabel('Node Centrality');title(plottitle);
    legend(legstr);
    mkNewDir(fullfile(outputpth,condition));
    mysave(gcf, fullfile(outputpth,condition,strcat(name,'3conditions', strsuffix)), 'all');
end



