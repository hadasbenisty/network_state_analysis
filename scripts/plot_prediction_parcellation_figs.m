addpath(genpath('../utils'));
addpath(genpath('../functions/'));
addpath(genpath('../meta_data_processing/'));
addpath(genpath('../svm'));
animal='xx';
statename =  'high_pup_l';
isloose = true;
loosestr = 'loose';
inputfolder = ['X:\Lav\ProcessingDirectory_Oct2020\' animal '\'];
outputfiggolder = 'C:\Users\Hadas Ben Esti\Dropbox (Personal)\cosyne2021\parcellation\figs\xx_loc\';
parcelsallen=load('X:\Hadas\Meso-imaging\Antara\preprocessing\parcells_updated121519.mat');
[parcels_names, ~, finalindex] = get_allen_meta_parcels;


%% Gal
load(strcat(inputfolder,'behavior_prediction',statename, loosestr ,'.mat'),'slopeData_gal',...
    'accuracy_vec_gal','accuracy_mat_gal');
acc_per_parcel_gal = scores_to_heatmap_gal(mean(accuracy_mat_gal), animal);
% acc_per_parcel_gal(isnan(acc_per_parcel_gal))=0.65;
figure;
% acc_per_parcel_gal(isnan(acc_per_parcel_gal))=-0.1;
% imagesc(acc_per_parcel_gal,[0.5 .8]);
pcolor([acc_per_parcel_gal nan(256,1); nan(1,256+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(ax, 'clim', [0.5 .8]);
colormap(jet);
axis off;c=colorbar;hold on;c.Limits=[0.5 0.8];
c.Ticks = [0.5 0.65 0.8];
c.FontSize=14;mysave(gcf, fullfile(outputfiggolder, 'gal_loc_acc_noblack'));

plot_parcellation_boundaries(parcelsallen.parcells_new.indicators(:,:,finalindex));
mysave(gcf, fullfile(outputfiggolder, 'gal_loc_acc'));
%% Allen
load(strcat(inputfolder,'behavior_prediction',statename, loosestr ,'.mat'),...
    'accuracy_mat');
acc_per_parcel_allen = scores_to_heatmap_allen(mean(accuracy_mat));
figure;pcolor([acc_per_parcel_allen nan(256,1); nan(1,256+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(ax, 'clim', [0.5 .8]);
;colormap(jet);

axis off;c=colorbar;
c.Limits=[0.5 0.8];
c.Ticks = [0.5 0.65 0.8];
c.FontSize=14;
hold on;mysave(gcf, fullfile(outputfiggolder, 'allen_loc_acc_noblack'));
plot_parcellation_boundaries(parcelsallen.parcells_new.indicators(:,:,finalindex));
mysave(gcf, fullfile(outputfiggolder, 'allen_loc_acc'));


Gv=[4 8 16];
for g=1:length(Gv)
    load(['X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_Slope_Amplitude\',animal,'\',animal,'trials_3states_grid' num2str(Gv(g)) loosestr '.mat'],...
        'parinds');
    [parcels_names_grid, parcels_region_labels_grid, final_index_grid, region_lut,...
        grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(parinds, Gv(g));
    
    load(strcat(inputfolder,'behavior_prediction',statename, '_grid', num2str(Gv(g)), loosestr,'.mat'), 'accuracy_mat');
    
    A = scores_to_heatmap_grid(accuracy_mat, animal, Gv(g), loosestr);
    figure;pcolor([A nan(256,1); nan(1,256+1)]);
    shading flat;
    set(gca, 'ydir', 'reverse');
    
    ;colormap(jet);
    
    
    
    axis off;c=colorbar;
    hold on;axis off;
    c.Limits=[0.5 0.8];
    c.Ticks = [0.5 0.65 0.8];
    c.FontSize=14;
    hold on;mysave(gcf, fullfile(outputfiggolder, ['grid_' num2str(Gv(g)) '_loc_acc_noblack']));
    
    plot_parcellation_boundaries(parcelsallen.parcells_new.indicators(:,:,finalindex));
    mysave(gcf, fullfile(outputfiggolder, ['grid_' num2str(Gv(g)) '_loc_acc']));
end










