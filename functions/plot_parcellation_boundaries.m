function plot_parcellation_boundaries(parcellsIndicators)
hold all;
for k=1:size(parcellsIndicators,3)
    E=parcellsIndicators(:,:,k);
    [B,~] = bwboundaries(E);hold on;
    for l = 1:length(B)
        boundary = B{l};
        plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1);
    end
end
end