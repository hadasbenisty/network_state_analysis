function W = process_sim(W1, th)

W = threshold_cor_matrix(W1, th);




% W=W-eye(size(W));
% W=W+1;
% W2=1-(W1+eye(size(W1)));
% W = gsp_learn_graph_log_degrees(W2, th, .1);
% W(W<1e-2) = 0;