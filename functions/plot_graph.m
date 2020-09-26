function plot_graph(G, cent, names, numeric_labels)

try
    figure;subplot(3,3,1);p=plot(G);
    p.NodeCData = numeric_labels;
    p.MarkerSize=4;
    for k=1:length(names)
        subplot(3,3,k+1);p=plot(G);
        p.MarkerSize=4;
        p.NodeCData = cent.(names{k});title(names{k});
    end
end
end