function [data_filt, W] = local_non_localmean(maskinds, R,C,data, patchsize, th, dispstr, timeinds)
if ~exist('timeinds','var')
    timeinds = 1:size(data,2);
end
data_filt = zeros(size(data));
[Y,X]=meshgrid(1:R,1:C);
X=X(:);X=X(maskinds);
Y=Y(:);Y=Y(maskinds);
% params = SetGenericDimsQuestParams(2, false);
% params.init_aff{1}.metric='cosine_similarityExp';
% params.init_aff{2}.metric='cosine_similarityExp';
% params.tree{1}.eigs_num=100;
% figure;subplot(2,2,1);
% P=zeros(R,C);
% P(maskinds) = data(:,1);
% imagesc(P);
figure;
W = spalloc(length(X),length(X),patchsize*patchsize*length(X));
tic;
for pii=1:length(X)
  
    x = X(pii);y = Y(pii);
    nn = find(abs(X-x) < patchsize & abs(Y-y) < patchsize);
    %     nn = setdiff(nn, pi);
    Cmat=corr(data(nn,timeinds).');
    w = Cmat(nn==pii,:);
    w(w < th) = 0;
    w=w/sum(w);
    data_filt(pii,:) = sum(bsxfun(@times, data(nn, :), w'),1);
    W(pii, nn) = w;
    if mod(pii,1000)==0
        disp(['Processed ' dispstr ' data: ' num2str(100*pii/length(X)) '%']);  
        toc;
        tic;
        P=zeros(R,C);
        P(maskinds(1:pii)) = data_filt(1:pii,1000);
        imagesc(P);drawnow;
      
    end
end
