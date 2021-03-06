function ii = discard_inds
load brain_mask.mat;
files = dir(['X:\Hadas\Meso-imaging\lan\' 'xx' 'psych\spt\' 'xx' '*_grid4.mat']);
load(fullfile(files(1).folder, files(1).name), 'par_inds');
[parcels_names_grid, parcels_region_labels_grid, final_index_grid, region_lut, grid_map_final_index, labelsbyallen] = getAllenClusteringLabelsGrid(par_inds, 4);
par_inds=par_inds(:, final_index_grid);
ii=[];
for k=1:size(par_inds,2)
x = par_inds(1,k)-4:par_inds(1,k)+4;
y = par_inds(2,k)-4:par_inds(2,k)+4;
if all(all(brain_mask(x,y)>0))
    ii(end+1) = k;
end
end
ii=1:921;
return;
disc = [9
    10
    11
    24
    25
    26
    41
    42
    43
    59
    62
    63
    83
    84
    85
   105
   106
   107
   130
   131
   154
   155
   178
   179
   203
   204
   228
   229
   254
   277
   278
   279
   302
   303
   304
   328
   329
   351
   352
   353
   354
   378
   379
   404
   454
   479
   530
   691
   712
   725
   732
   745
   746
   752
   771
   785
   786
   787
   788
   789
   790
   802
   803
   804
   805
   806
   807
   808
   809
   819
   821
   823
   824
   825
   826
   832
   834
   835
   836
   837
   838
   839
   840
   841
   842
   843
   844
   845
   846
   848
   850
   852
   853
   854
   855
   856
   857
   858
   859
   860
   861
   862
   863
   864
   866
   867
   868
   869
   870
   871
   872
   873
   874
   875
   876
   877
   878
   879
   880
   881
   882
   883
   884
   885
   886
   887
   888
   889
   890
   891
   892
   894
   895
   896
   897
   898
   899
   901
   902
   904
   905
   906
   907
   908
   909
   910
   911
   912
   913
   914
   915
   917
   918
   919
   920
   921];

ii=setdiff(1:921,disc);