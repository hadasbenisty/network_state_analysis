function [indic, cent, G, c] = graph_analysis(W, th, str_labels)
if ~exist('th','var')
    th = 0.01;
end
W(W<th)=0;
isdirected = ~issymmetric(W);
if isdirected
    G=digraph(W, str_labels);
    c = {'indegree','outdegree','incloseness','outcloseness','betweenness','pagerank','hubs','authorities'};
    
else
    G = graph(W);
    c = {'degree','closeness','betweenness','pagerank', 'eigenvector'};
end
indic = zeros(length(c), size(W, 1));
for k=1:length(c)
    cent.(c{k}) = centrality(G, c{k});
    indic(k,(cent.(c{k})> median(cent.(c{k}))))=1;
end