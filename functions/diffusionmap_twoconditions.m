function diffusionmap_twoconditions(W1,W2,centrality1,centrality2,name1,name2,Session,centname)
for i=1:2
    if i==1
        W=W1;
        centrality=centrality1;
        name=name1;
    elseif i==2
        W=W2;
        centrality=centrality2;
        name=name2;
    end
    %threshold the weight matrix
    params.knn =5;
    nn_dist = sort(abs(W.'),'descend').';
    v = nn_dist(:, 1:params.knn);
    th= median(v(:));
    W(abs(W)<th)=0;
    K = 1-(W+W')/2;
    K = K.*(1-eye(size(K)));
    %get diffusion map
    [diffusion_map, Lambda, Psi, Ms, Phi, K_rw] = calcDiffusionMap(K,[]);
    %make labels for plot
    toremove=setdiff(1:56,[21:26 53:56]);isleftlabel=2:2:56;
    finalindex=intersect(isleftlabel,toremove);
    parcels_region_labels_bilateral=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 0 0 0 0 3 3 4 4 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 0 0 0 0];
    parcels_region_labels=parcels_region_labels_bilateral(finalindex);
    %Plot 1
    figure;plotEmbeddingWithColors(diffusion_map', parcels_region_labels,name);
    title(strcat(name,Session));
    set(gcf,'renderer','Painters');
    mysave(gcf,fullfile('X:\Hadas\Meso-imaging\lan\results\ProcessingDirectory\allen_directed_spon_network\',Session,strcat('\','regionlabel',name)),'all');%
    clearvars -except W1 W2 centrality1 centrality2 name1 name2 Session centname min_c max_c
end
end

