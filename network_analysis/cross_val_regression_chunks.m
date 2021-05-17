function [rsquare, yfit] = cross_val_regression_chunks(kfold, X, y)
[n,nt,nT] = size(X);

if size(y, 1) ~= nt || size(y,2) ~= nT
    error('dims of X and y are incompatible');
end

if nT > 1
    nt = size(X,2);
    X = reshape(X, n, nt*nT);
end
y=y(:);
yfit = nan(round(kfold*nt), 1);
for ii=1:kfold
    teinds = 1+(ii-1)*nt:ii*nt;
    trinds = setdiff(1:kfold*nt,teinds);
    if any([trinds teinds] > size(X,2))
        break;
    end
    mdl = fitlm(X(:,trinds)',y(trinds));
    yfit(teinds) = predict(mdl,X(:,teinds)');
end
yfit=yfit(~isnan(yfit));
yfit=yfit(:);

SSE = mean((yfit-y(1:length(yfit))).^2);
TSS = var(y(1:length(yfit)));
rsquare = 1 - SSE/TSS;

yfit=reshape(yfit, nt, nT);