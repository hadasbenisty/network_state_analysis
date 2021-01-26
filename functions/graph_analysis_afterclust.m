function [cent_weighted, cent_notweighted, G, c] = graph_analysis_afterclust(W1, str_labels, communitylabels,processWfun, th, issim)
addpath(genpath('../centrality_measures/'));
addpath(genpath('..\gspbox'));
if ~exist('processWfun','var')
processWfun = @process_sim;
% processWfun = @id;
end
if ~exist('th','var')
th = 7;
end
W = feval(processWfun, W1, th);


if ~exist('issim','var')
    issim=1;
end
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
isdirected = ~issymmetric(W) || ~issim;

if isdirected
    G=digraph(double(W), str_labels);
    c = {'indegree',             'outdegree',        'incloseness','outcloseness',...
        'betweenness','pagerank','hubs','authorities'};
    wstr = {'Importance'   'Importance'               'Cost'      'Cost'    'Cost'      'Importance'   'Importance' 'Importance'};
    wnums = {abs(G.Edges.Weight) abs(G.Edges.Weight) 1./abs(G.Edges.Weight) 1./abs(G.Edges.Weight)...
         1./abs(G.Edges.Weight) abs(G.Edges.Weight) abs(G.Edges.Weight) abs(G.Edges.Weight)};
else
    G = graph(W);
    c = {'degree',         'closeness','betweenness','pagerank', 'eigenvector'};
    wstr = {'Importance'   'Cost'      'Cost'         'Importance'   'Importance'};
    wnums = {abs(G.Edges.Weight) 1./abs(G.Edges.Weight) 1./abs(G.Edges.Weight)...
        abs(G.Edges.Weight) abs(G.Edges.Weight)};
end

for k=1:length(c)    
    cent_weighted.(c{k}) = centrality(G, c{k}, wstr{k}, wnums{k}+eps);    
    cent_notweighted.(c{k}) = centrality(G, c{k});
end
if exist('communitylabels','var') && ~isempty(communitylabels)
cent_notweighted.participation=participation_coef((W~=0),communitylabels,0);
cent_notweighted.community = communitylabels;

cent_weighted.participation=participation_coef(W,communitylabels,0);
cent_weighted.community = communitylabels;
end
configParams.maxInd=10;
[diffusion_map, Lambda] = calcDiffusionMap(W,configParams);
cent_weighted.diffmap = diffusion_map(1,:)';
cent_weighted.second_eigval = 1-Lambda(2);
[diffusion_map, Lambda] = calcDiffusionMap(W>0,configParams);
cent_notweighted.diffmap = diffusion_map(1,:)';
cent_notweighted.second_eigval = 1-Lambda(2);

end