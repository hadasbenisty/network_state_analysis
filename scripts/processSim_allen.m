function aff_mat=processSim_allen(W,th)

params.knn =th;

nn_dist = sort(abs(W.'),'descend').';
v = nn_dist(:, 2:params.knn);
th= median(v(:));
aff_mat=abs(W);
aff_mat(abs(W)<th)=0;
% aff_mat=aff_mat.*(1-eye(size(W)));