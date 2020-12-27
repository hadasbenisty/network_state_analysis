function aff_mat = threshold_cor_matrix(W,k)

params.knn =k;
% if size(W,1)>100
% %     params.knn =100;
%     nn_dist = sort(abs(W.'),'descend').';
%     v = nn_dist(:, 1:params.knn);
%     th= 3*median(v(:));
% else
    
    nn_dist = sort(abs(W.'),'descend').';
    v = nn_dist(:, 1:params.knn);
    th= median(v(:));
% end
aff_mat=abs(W);
aff_mat(abs(W)<th)=0;
% params.knn =round(size(W,1)/3);

% params.knn =5;
% nn_dist = sort(abs(W.'),'descend').';
% v = nn_dist(:, 1:params.knn);
% th= median(v(:));
% aff_mat=abs(W);
% aff_mat(abs(W)<th)=0;
% nn_dist = sort(1-W.').';
% params.knn = 5;
% sigma = 0.1 * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
% params.knn=10;
% aff_mat = exp(W/(sigma));
% nn_dist = sort(aff_mat.','ascend').';
% v = nn_dist(:, 2:params.knn+1);
% th= 10*median(v(:));
% 
% aff_mat(abs(aff_mat)<th)=0;






% params.knn =round(size(W,1)/3);

% params.knn =10;
% % W=abs(W);
% nn_dist = sort(abs(W).','ascend').';
% % v = nn_dist(:, 2:params.knn+1);
% % th= 10*median(v(:));
% th= median(nn_dist(:));
% aff_mat=W;
% aff_mat(abs(W)<th)=0;

% 
% params.knn =10;
% nn_dist = sort(abs(W).','ascend').';
% v = nn_dist(:, 2:params.knn+1);
% th= 10*median(v(:));
% aff_mat=(exp(abs(W)/th));

% params.knn =20;
% nn_dist = sort(aff_mat.','ascend').';
% v = nn_dist(:, 2:params.knn+1);
% th= median(nn_dist(:));
% aff_mat(abs(aff_mat)<th)=0;
