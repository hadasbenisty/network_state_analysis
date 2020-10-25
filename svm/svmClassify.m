%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% svm prediction - training and testing based on libsvm
% input:
%   X               - d x n matrix of features where n is the amount of
%                     vectors and d is their dimension
%   Y               - n labels
%   foldsNum        - folds for cross validation
%   tonorm          - whether or not to normalize data set
%   kernelType      -
%                       0 -- linear: u'*v
%                       1 -- polynomial: (gamma*u'*v + coef0)^degree
%                       2 -- radial basis function: exp(-gamma*|u-v|^2)
%                       3 -- sigmoid: tanh(gamma*u'*v + coef0)
%   bootstrapmethod - a string indicating bootstapping method - 'none',
%                     'downsample' or 'upsample'
% output:
%   accuracy_vec    - a vector consisting the accuracy rate for each fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accuracy_vec, cvinds, fa_vec, md_vec, accuracy_vecSh, fa_vecSh, md_vecSh]  = svmClassify(X, Y, foldsNum, kernelType, tonorm, bootstrapmethod, upsamplesz, voting, doshuffle, REPS)
if ~exist('voting','var')
    voting = false;
end
if ~exist('doshuffle','var')
    doshuffle = false;
end
cvinds=[]; fa_vec=[]; md_vec=[];
% handle input vars
Y = Y(:);
if size(X,1) ~= length(Y)
    error('Dims are inconsistent');
end

if tonorm
    Xnorm = (X - min(X(:)))/(max(X(:))-min(X(:)));
else
    Xnorm=X;
end
Xnorm = Xnorm.';
log2c = -6:10;log2g = -6:4;


kernelStr = [' -t ' num2str(kernelType)];
if kernelType == 0
    log2g = [];
end
% to use for bootstrapping with upsample (if chosen)
if ~exist('upsamplesz', 'var')
    upsamplesz = 100;
end

% output var
accuracy_vec = zeros(foldsNum, 1);

clusters = unique(Y);
% evaluate weight of each cluster - used only without bootstrapping
if strcmp(bootstrapmethod, 'none')
    wstr=' ';
    for w_i=1:length(clusters)
        w = sum(Y==clusters(w_i))/length(Y);
        wstr = [wstr ' -w' num2str(clusters(w_i)) ' ' num2str(w)]; %#ok<AGROW>
    end
else
    wstr = '';
end
wstr = [wstr ' '];

% draw cross validation partition
clustersSz = hist(Y, unique(Y));
if any(clustersSz < 2*foldsNum) || length(clustersSz) == 1
    warning('Cannot cross validate, too few examples');
    accuracy_vec = [];return;
end
for ci = 1:length(clusters)
    clusterinds = find(Y == clusters(ci));
    cvindspercluster{ci} = crossvalind('Kfold', length(clusterinds), foldsNum);
    cvinds(clusterinds) = cvindspercluster{ci};
end

% cross-validation loop for training and testing
for foldi = 1:foldsNum
    % divide train and test sets
    testinds = find(cvinds == foldi);
    traininds = setdiff(1:length(Y), testinds);
    Xtrain = Xnorm(:, traininds);
    Ytrain = Y(traininds);
    Xtest  = Xnorm(:, testinds);
    Ytest   = Y(testinds);
    if length(unique(Ytest)) ~= length(unique(Ytrain))
        error('At least one cluster is too small');
    end
    
    switch bootstrapmethod
        case 'none'
            XcurrTr = Xtrain;
            YcurrTr = Ytrain;
            XcurrTe = Xtest;
            YcurrTe = Ytest;
        case 'downsample'
            % find the smallest cluster
            smallestClassSzTraining = min(hist(Ytrain, unique(Ytrain)));
            smallestClassSzTesting  = min(hist(Ytest, unique(Ytest)));
            XcTr = cell(length(clusters), 1);
            XcTe = cell(length(clusters), 1);
            SzTr = zeros(length(clusters), 1);
            SzTe = zeros(length(clusters), 1);
            XcurrTr=[];XcurrTe=[];
            for ci = 1:length(clusters)
                XcTr{ci} = Xtrain(:, Ytrain == clusters(ci));
                SzTr(ci) = sum(Ytrain == clusters(ci));
                XcTe{ci} = Xtest(:, Ytest == clusters(ci));
                SzTe(ci) = sum(Ytest == clusters(ci));
                randinds = randperm(SzTr(ci));
                XcurrTr = cat(2, XcurrTr, XcTr{ci}(:, randinds(1:smallestClassSzTraining)));
                randinds = randperm(SzTe(ci));
                XcurrTe = cat(2, XcurrTe, XcTe{ci}(:, randinds(1:smallestClassSzTesting)));
            end
            YcurrTr = kron(clusters,ones(smallestClassSzTraining, 1));
            YcurrTe = kron(clusters,ones(smallestClassSzTesting, 1));
        case 'upsample'
            smallestClassSzTesting  = min(hist(Ytest, unique(Ytest)));
            XcTr = cell(length(clusters), 1);
            XcTe = cell(length(clusters), 1);
            SzTr = zeros(length(clusters), 1);
            SzTe = zeros(length(clusters), 1);
            XcurrTr=[];XcurrTe=[];
            for ci = 1:length(clusters)
                XcTr{ci} = Xtrain(:, Ytrain == clusters(ci));
                SzTr(ci) = sum(Ytrain == clusters(ci));
                XcTe{ci} = Xtest(:, Ytest == clusters(ci));
                SzTe(ci) = sum(Ytest == clusters(ci));
                indsTr = randsample(SzTr(ci), upsamplesz, true);
                XcurrTr = cat(2, XcurrTr, XcTr{ci}(:, indsTr));
                randinds = randperm(SzTe(ci));
                XcurrTe = cat(2, XcurrTe, XcTe{ci}(:, randinds(1:smallestClassSzTesting)));
            end
            YcurrTr = kron(clusters,ones(upsamplesz, 1));
            YcurrTe = kron(clusters,ones(smallestClassSzTesting, 1));
    end
    
    if voting
        [accuracy_vec(foldi), ~, ~, fa_vec(:, foldi), md_vec(:, foldi)] = trainPredictGetAccVoting(XcurrTr.', YcurrTr, log2c, log2g, foldsNum, kernelStr, wstr, XcurrTe.', YcurrTe);
    else
    [accuracy_vec(foldi), ~, ~, fa_vec(:, foldi), md_vec(:, foldi)] = trainPredictGetAcc(XcurrTr.', YcurrTr, log2c, log2g, foldsNum, kernelStr, wstr, XcurrTe.', YcurrTe);
    if doshuffle
        for ri = 1:REPS
       XcurrTr1 = shuffleFeatures(XcurrTr, YcurrTr);
       XcurrTe1 = shuffleFeatures(XcurrTe, YcurrTe);
       [accuracy_vec1(foldi, ri), ~, ~, fa_vec1(:, foldi, ri), md_vec1(:, foldi, ri)] = trainPredictGetAcc(XcurrTr1.', YcurrTr, log2c, log2g, foldsNum, kernelStr, wstr, XcurrTe1.', YcurrTe);
        end
    end
    
    end
    
end
if doshuffle
accuracy_vecSh = mean(accuracy_vec1,2);
fa_vecSh = mean(fa_vec1,3);
md_vecSh = mean(md_vec1,3);
else
    accuracy_vecSh=[];
    fa_vecSh=[];
md_vecSh=[];
end
end

function [acc, SVMModel, score, fa, md] = trainPredictGetAcc(Xtrain, Ytrain, log2c, log2g, foldsNum, kernelStr, wstr, Xtest, Ytest)

% use training data to determine parameters
[~, ~, cvparamsStr] = cvsvmclassification(Xtrain, Ytrain, 2.^log2c, 2.^log2g, foldsNum);
% train model with best parameters
SVMModel = svmtrain(Ytrain, Xtrain, [kernelStr ' -q ' cvparamsStr wstr]); %#ok<*SVMTRAIN>
% predict
[predictions, ~, score] = svmpredict(Ytest, sparse(Xtest), SVMModel);
% evaluate accuracy
acc = sum(predictions==Ytest)/length(Ytest);
classes = unique(Ytest);
for ci = 1:length(classes)
    fa(ci) = sum((predictions == classes(ci)) & (Ytest ~= classes(ci)))/length(Ytest);
    md(ci) = sum((predictions ~= classes(ci)) & (Ytest == classes(ci)))/length(Ytest);
end


end


function [acc, SVMModel, score, fa, md] = trainPredictGetAccVoting(Xtrain, Ytrain, log2c, log2g, foldsNum, kernelStr, wstr, Xtest, Ytest)

% use training data to determine parameters
for voter_i = 1:size(Xtrain,2)
[~, ~, cvparamsStr] = cvsvmclassification(Xtrain(:,voter_i), Ytrain, 2.^log2c, 2.^log2g, foldsNum);
% train model with best parameters
SVMModel = svmtrain(Ytrain, Xtrain(:,voter_i), [kernelStr  cvparamsStr wstr ' -q ']); %#ok<*SVMTRAIN>
% predict
[predictions(:,voter_i),~, score] = svmpredict(Ytest, sparse(Xtest(:,voter_i)), SVMModel);
end
classes = unique(Ytest);
for p_i = 1:size(predictions,1)
b = hist(predictions(p_i,:), unique(Ytest));
[~,m] = max(b);
finalPred(p_i) = classes(m);
end
finalPred=finalPred(:);
% evaluate accuracy
acc = sum(finalPred==Ytest)/length(Ytest);
classes = unique(Ytest);
for ci = 1:length(classes)
    fa(ci) = sum((finalPred == classes(ci)) & (Ytest ~= classes(ci)))/length(Ytest);
    md(ci) = sum((finalPred ~= classes(ci)) & (Ytest == classes(ci)))/length(Ytest);
end


end





