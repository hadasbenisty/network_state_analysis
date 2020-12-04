function W = threshold_cor_matrix(W)
params.knn =round(size(W,1)/3);

% params.knn =5;
nn_dist = sort(abs(W.'),'descend').';
v = nn_dist(:, 1:params.knn);
th= median(v(:));
W(abs(W)<th)=0;


