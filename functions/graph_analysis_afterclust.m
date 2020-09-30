function [indic, cent, G, c] = graph_analysis_afterclust(W, str_labels, isweighted)
params.knn =5;
%step 1
%W=1-W./(mean(W(:)).^2);
%W=1-W.^2;
%step 2
%nn_dist = sort(W.').';

nn_dist = sort(abs(W.'),'descend').';
v = nn_dist(:, 1:params.knn);
th= median(v(:));%.10l*mean(v(v~=0));
W(abs(W)<th)=0;

%% adaptive
% Wthresholded=zeros(23,23);
% Windx=zeros(23,23);
% for j=1:23
%     tempvec=W(j,:);
%     tempvecindx=tempvec<th(j); % if below average correlation, store as 1.
%     tempvec(tempvec<th(j))=0;
%     Wthresholded(j,:)=tempvec;
%     Windx(j,:)=1-tempvecindx; % if below average correlation, do not index (reverse of logical statement from earlier
% end
%W=Wthresholded;

%step 3
% expvalue=bsxfun(@rdivide, -Wthresholded.^2, th.^2.');
% W=exp(expvalue);
%W(Windx==0)=0;
%W = bsxfun(@rdivide, W, sum(W,2)+eps);
%W(eye(size(W)) == 1) = 0;
%step 4
isdirected = ~issymmetric(W);

if isdirected
    G=digraph(W, str_labels);
    c = {'indegree','outdegree','incloseness','outcloseness','betweenness','pagerank','hubs','authorities'};
else
    G = graph(W);
    c = {'degree','closeness','betweenness','pagerank', 'eigenvector'};
end
indic = zeros(length(c), size(W, 1));


if isweighted
    for k=1:length(c)
        switch c{k}
            case {'degree', 'eigenvector' , 'pagerank'}
                cent.(c{k}) = centrality(G, c{k}, 'Importance',exp(abs(G.Edges.Weight)));
            case { 'closeness' 'betweenness'}
                cent.(c{k}) = centrality(G, c{k}, 'Cost',exp(-G.Edges.Weight));
        end
        indic(k,(cent.(c{k})> median(cent.(c{k}))))=1;
    end
else
    for k=1:length(c)
        cent.(c{k}) = centrality(G, c{k});
        indic(k,(cent.(c{k})> median(cent.(c{k}))))=1;
    end
end