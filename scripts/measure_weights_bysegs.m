%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measures strength of connectivity of time traces.
% Input -
%   * data - p time traces over T time frames
%   * method - how to measure weights
%   * params - a struct of internal parameters for each method
%
% Output - weights matrix
%
% Written by: Hadas Benisty, July 8, 2020
% hadas.benisty@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [W, Wmat, well_modeled_nodes] = measure_weights_bysegs(t, data, method, params)
t=t(1:size(data,2));
if isempty(data)
    W=[];
    return;
end
well_modeled_nodes = [];
if ~exist('method', 'var')
    method = 'corr';
    params.is_symmetric = true;
    params.zero_diag = true;
    params.how_to_normalize = true;
    params.how_to_normalize = 'rows';
end
if ~exist('param', 'var')
    switch method
        case 'relative_modeling_contribution'
            params.is_symmetric = false;
            params.zero_diag = true;
            params.how_to_normalize = true;
            params.how_to_normalize = 'rows';
            
            params.folds_num = 5;
            params.model = 'lasso';
            params.modeled_energy_th = 0.1;
        case {'partial_corr','partial_corr_mean_pop'}
            params.is_symmetric = true;
            params.zero_diag = true;
            params.how_to_normalize = false;
            params.how_to_normalize = 'rows';
        case {'pearson_corr','cov'}
            params.is_symmetric = true;
            params.zero_diag = true;
            params.how_to_normalize = false;
            params.how_to_normalize = 'rows';
        case {'L2'}
            params.is_symmetric = true;
            params.zero_diag = true;
            params.how_to_normalize = false;
            params.how_to_normalize = 'rows';
    end
end
Td = (find(diff(t)>1.1*median(diff(t))));
Td=[Td; length(t)];
st=1;
for i = 1:length(Td)
    y = data(:, st:Td(i));
switch method
    case 'pearson_corr'
        W = corr(y');
    case 'cov'
        y=bsxfun(@minus, X, mean(y,2));
        inner_products     = y * y.';
        [ij_norm, ji_norm] = meshgrid( diag(inner_products) );
        W            = inner_products ./ (sqrt(ij_norm .* ji_norm)+eps);
           

case 'partial_corr_mean_pop'
        %W = abs(corr(data'));
        pairs = nchoosek(1:size(y, 1), 2);
        rho=[];
        for i=1:size(pairs,1)
            X=y(pairs(i,1),:).';
            Y=y(pairs(i,2),:).';
            Z=y(setdiff(1:size(y,1),[pairs(i,1),pairs(i,2)]),:).'; %for partial correlation, all parcels outside of the pair to control for
            rho(i)=partialcorr(X,Y,nanmean(Z,2));
            clear vars X Y Z
        end
        %        C=corr(xet');
        %         C = C.*(C>th);
        W=squareform(rho);
        W(eye(size(W)) == 1) = 1;

    case 'partial_corr'
        %W = abs(corr(data'));
        pairs = nchoosek(1:size(data, 1), 2);
        rho=[];
        for i=1:size(pairs,1)
            X=data(pairs(i,1),:).';
            Y=data(pairs(i,2),:).';
            Z=y(setdiff(1:size(data,1),[pairs(i,1),pairs(i,2)]),:).'; %for partial correlation, all parcels outside of the pair to control for
            rho(i)=partialcorr(X,Y,Z);
            clear vars X Y Z
        end
        %        C=corr(xet');
        %         C = C.*(C>th);
        W=squareform(rho);
        W(eye(size(W)) == 1) = 1;
    case 'relative_modeling_contribution'
        R2_full = train_full_model(y, params);
        well_modeled_nodes = find(R2_full > params.modeled_energy_th);
        R2_partial = train_partial_model(data, params);
        W = calc_contribution(R2_full, R2_partial, params.modeled_energy_th);
    case 'L2'
        dParams.kNN = 20;
        dParams.self_tune = 0;
        dParams.verbose = false;
        dParams.isuper = 1; %whether to perform supercharging in rann32 (0 = no)
        
        W = calcAffinityMat(y.', dParams);
end

if params.is_symmetric
    W=(W+W')/2;
end

if params.zero_diag
    W(eye(size(W)) == 1) = 0;
end
%
% if params.how_to_normalize
%     switch params.how_to_normalize
%         case 'rows'
%             W = bsxfun(@rdivide, W, sum(W,2)+eps);
%         case 'cols'
%             W = bsxfun(@rdivide, W, sum(W,1));
%     end
Wmat(:,:,i) = W;
end
W = nanmean(Wmat,3);
if any(isnan(W(:)))
    W=[];
end

end



