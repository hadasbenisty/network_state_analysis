function graph_overlay_allen_notpaired(outputpth, cond1_cat,cond2_cat,condition,name,plottitle,parcels_names,num)
    %calculate average and error of condition 1
    condition1_meanvec_parcel=mean(cond1_cat,2);
    condition1_sevec_parcel=std(cond1_cat,0,2)./sqrt(num-1);
    %calculate average and error of condition 2
    condition2_meanvec_parcel=mean(cond2_cat,2);
    condition2_sevec_parcel=std(cond2_cat,0,2)./sqrt(num-1);
    
    %get colors
    isleftlabel=2:2:56;
    toremove=setdiff(1:56,[21:26 53:56]);
    finalindex=intersect(isleftlabel,toremove);
    parcels_region_labels=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
    parcels_region_labels=parcels_region_labels(finalindex);
    RegionColors=[0 0.4470 0.7410;1 1 0.1294;0.4706 0.6706 0.3020;0.8 1 0;0 0.6 0.9;0.9294 0.6902 0.1294;0.9294 0.4000 0.1294;0.9294 0.1020 0.1294];
    %plot bar in loop to custom color
    figure;
    set(gcf,'renderer','Painters')
    for i=1:23
        indx=parcels_region_labels(i);
        h=errorbar(i,condition1_meanvec_parcel(i),condition1_sevec_parcel(i),'color',RegionColors(indx,:),'LineStyle','none','LineWidth',1.5);ylabel('Node Centrality');
        title(underline2space(plottitle));        
        alpha = 1;
        set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
        hold on
    end
        hold on    
    for i=1:23
        indx=parcels_region_labels(i);
        h=errorbar(i,condition2_meanvec_parcel(i),condition2_sevec_parcel(i),'color',RegionColors(indx,:),'LineStyle','none','LineWidth',1.5);ylabel('Node Centrality');
        title(underline2space(plottitle));        
        alpha = 0.3;
        set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
        hold on
    end
    set(gcf, 'Position',  [150,150, 1000,500]);
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    
    mkNewDir(condition);
    mysave(gcf, fullfile('X:\Lav\ProcessingDirectory\parcor_undirected\',condition,strcat(name,'regioncolored')), 'all');
    %plot the same thing again but color coded for condition not region
    figure;
    set(gcf,'renderer','Painters')
    errorbar(condition1_meanvec_parcel,condition1_sevec_parcel,'color',[0 0 0.8008],'LineStyle','none','LineWidth',1.5);ylabel('Node Centrality');title(underline2space(plottitle));
    hold on;
    errorbar(condition2_meanvec_parcel,condition2_sevec_parcel,'color',[0.9961 0.2695 0],'LineStyle','none','LineWidth',1.5);ylabel('Node Centrality');title(underline2space(plottitle));
    set(gcf, 'Position',  [150,150, 1000,500]);
    set(gca,'xtick',1:23)
    set(gca,'xticklabel',parcels_names)
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',15)
    set(gca,'XTickLabelRotation',45);
    mkNewDir(fullfile(outputpth,condition));
    mysave(gcf, fullfile(outputpth,condition,strcat(name,'conditioncolored')), 'all');
end

