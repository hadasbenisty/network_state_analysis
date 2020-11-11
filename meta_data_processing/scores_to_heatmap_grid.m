function A = scores_to_heatmap_grid(Xvec, animal, G, loosestr)

load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states_grid' num2str(G) loosestr '.mat'],...
    'parinds');
[~, ~, ~, ~,grid_map_final_index] = getAllenClusteringLabelsGrid(parinds, G);

c=unique(grid_map_final_index);
A=nan(size(grid_map_final_index));
v=mean(Xvec);
for i=2:length(c)
    A(grid_map_final_index==c(i)) = v(i-1);
    
end





