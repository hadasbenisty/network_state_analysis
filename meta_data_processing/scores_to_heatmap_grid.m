function A = scores_to_heatmap_grid(Xvec, animal, G)
files = dir(['X:\Hadas\Meso-imaging\lan\' animal 'psych\spt\' animal '*_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');

[~, ~, ~, ~,grid_map_final_index] = getAllenClusteringLabelsGrid(par_inds, G);

c=unique(grid_map_final_index);
A=nan(size(grid_map_final_index));

for i=2:length(c)
    A(grid_map_final_index==c(i)) = Xvec(i-1);
    
end





