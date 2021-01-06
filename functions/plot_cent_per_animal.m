function plot_cent_per_animal(heatmaps)
L = quantile(heatmaps(:),[.1 .9]);

figure;k=1;
set(gcf,'Position',[1          41        1920         963]);

for state_i = 1:size(heatmaps,4)
    
    for ai=1:size(heatmaps,3)
        subplot(size(heatmaps,4), size(heatmaps,3)+1, k);
%         x=heatmaps(:,:,ai,state_i);
%         L = quantile(x(:),[.1 .9]);

        
        plot_vals_heatmap(heatmaps(:,:,ai,state_i), 'Centrality',...
            [],  L(1)-eps, L(2)+eps, 1,colormap(redblue))
        title([' animal #' num2str(ai)  ]);k=k+1;
    end
    subplot(size(heatmaps,4), size(heatmaps,3)+1, k);
    
    
    plot_vals_heatmap(nanmean(heatmaps(:,:,:,state_i),3), 'Centrality',...
        [],   L(1)-eps, L(2)+eps, 1,colormap(redblue))
    title('Mean');k=k+1;
end
