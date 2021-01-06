function [ aff_mat, slices ] = CalcInitAff( data, params, dim )
% Calculates the affinity between slices by columns of the data as thresholded
% cosine similarity between columns, or as a Gaussian kernel of width
% eps*(median distance between 5 nearest neighbors of all points).
%
% Inputs:
%     data - M-by-N-by-nT matrix
%     params - struct array with user parameters
%
% Output:
%     aff_mat - N-by-N symmetric matrix of non-negative affinities
%--------------------------------------------------------------------------
dimLen = length(size(data));
% backwards compatibility
if ~isfield(params, 'RangeMinus1to1')
    params.RangeMinus1to1 = false;
end

data   = permute(data, [setdiff(1:dimLen, dim), dim]);
slices = reshape(data, [], size(data, dimLen));

if dimLen >= 3
    %     slices = zeros(size(data, 1)* size(data, 3), size(data, 2));
    %     for r=1:size(data, 2)
    %         currslice = data(:, r, :);
    %         slices(:, r) = currslice(:);
    %     end
    switch params.metric   
        case 'inputAff'
            aff_mat = params.input;
        case 'removeSideInfo'
             euc_dist = squareform(pdist(slices.'));
             euc_distSideInfo = squareform(pdist(params.sideinfo));
             euc_dist=euc_dist./(1+0.10*euc_distSideInfo);
            nn_dist = sort(euc_dist.').';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-euc_dist.^2/(2*sigma^2));
            
        case 'sideInfo'
            euc_dist = squareform(pdist(params.loc));
            nn_dist = sort(euc_dist.').';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-euc_dist.^2/(2*sigma^2));
           
            
         case 'covRienmanMean'
            for nr = 1:size(data,3)
                C(:,:, nr) = cov(data(:, :, nr));
            end
            
            for nr1 = 1:size(data,3)
                for nr2 = 1:size(data,3)
                    meanC = riemann_mean(C(:,:, [nr1, nr2]));
                    aff_mat(nr1, nr2) = distance_riemann(meanC, C(:,:, nr1))+distance_riemann(meanC, C(:,:, nr2));
%                     aff_mat(nr1, nr2) = RiemannianDist(meanC, C(:,:, nr2));
                end
            end   
        case 'covRienman'
            for nr = 1:size(data,3)
                C(:,:, nr) = cov(data(:, :, nr));
            end
            for nr1 = 1:size(data,3)
                for nr2 = 1:size(data,3)
                    aff_mat(nr1, nr2) = RiemannianDist(C(:,:, nr1), C(:,:, nr2));
                end
            end
            
        case 'cosine_similarityOnTrials'
            for T=1:size(data, 3)
                v = permute(data(:, :, T), [1 2 3]);
                inner_products = v.'*v;
                [ij_norm, ji_norm] = meshgrid(diag(inner_products));
                aff_matAll(:,:,T) = inner_products./sqrt(ij_norm.*ji_norm);
            end
            aff_mat = mean(aff_matAll, 3);
            aff_mat(aff_mat < params.thresh) = 0;
        case 'corr'
            aff_mat     = slices.' * slices;
            aff_mat(aff_mat < params.thresh) = 0;
        case 'cosine_similarity'
            inner_products     = slices.' * slices;
            [ij_norm, ji_norm] = meshgrid( diag(inner_products) );
            aff_mat            = inner_products ./ (sqrt(ij_norm .* ji_norm)+eps);
            aff_mat(aff_mat < params.thresh) = 0;
        case 'cosine_similarityExp'
            inner_products     = slices.' * slices;
            [ij_norm, ji_norm] = meshgrid( diag(inner_products) );
            aff_mat            = inner_products ./ (sqrt(ij_norm .* ji_norm)+eps);
%             aff_mat(aff_mat < params.thresh) = 0;
            sigma = params.eps * median(aff_mat(:));
            if sigma == 0
                sigma = params.eps * mean(aff_mat(:));
            end
            aff_mat = exp(aff_mat/(sigma));
            
        case 'emd'
            [L1, L2] = size(slices);
            for n1=1:L2
                for n2 = 1:L2
                    [~, aff_mat(n1, n2)] = emd((1:L1).', (1:L1).', slices(:,n1), slices(:,n2), @gdf);
                end
            end
        case 'L1'
            
            
            euc_dist = squareform(pdist(slices.', 'cityblock'));
            nn_dist = sort(euc_dist.').';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-euc_dist/(sigma));
        case 'dtw'
            
            for k1 = 1:size(slices,2)
                for          k2 = k1+1:size(slices,2)
                    dtw_dist(k1,k2) = dtw_c(slices(:, k1), slices(:, k2), params.maxDTWlen);
                end
            end
            for k1 = 1:size(slices,2)
                for          k2 = k1+1:size(slices,2)
                    dtw_dist(k2,k1) = dtw_dist(k1,k2);
                end
            end
            nn_dist = abs(sort(dtw_dist.')).';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-dtw_dist.^2/(2*sigma^2));
        case {'euc','L2'}
            
            
            euc_dist = squareform(pdist(slices.'));
            nn_dist = sort(euc_dist.').';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-euc_dist/sigma);
        case 'spikeTrainByVictor'
            for ind3 = 1:size(data, 2)
                disp(['spikeTrainByVictor: third dim:' num2str(ind3)]);
                [ aff_mat(:,:,ind3) ] = CalcInitAff2D( squeeze(data(:,ind3,:)), params );
            end
            aff_mat = mean(aff_mat, 3);
        case 'connectivity'
            aff_mat = GraphConnectivity( permute(data, [1 3 2]),  params);
        case {'weuc','wL2'}
            
            w = mean(data,3).';
            w=w(:)/sum(w(:));
            wslices = bsxfun(@times, slices, w);
            euc_dist = squareform(pdist(wslices.'));
            nn_dist = sort(euc_dist.').';
            params.knn = min(params.knn, size(nn_dist, 2));
            sigma = params.eps * median(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            if sigma == 0
                sigma = params.eps * mean(reshape(nn_dist(:, 1:params.knn), size(nn_dist,1)*params.knn,1));
            end
            aff_mat = exp(-euc_dist/sigma);
        otherwise
            error('params.metric is not well define, please fix it');
            
    end
    
    if params.RangeMinus1to1
        aff_mat = 2*(aff_mat-min(aff_mat(:)))./(max(aff_mat(:)-min(aff_mat(:))))-1;
    end
else
    aff_mat  = CalcInitAff2D( data, params );
    
    % figure, imagesc(aff_mat), colormap gray, colorbar, axis on, hold on
    % if params.on_rows,
    %     title('Initial Row Affinity'), hold off
    % else
    %     title('Initial Col Affinity'), hold off
    % end
    
end
if isfield(params, 'mask')
    aff_mat=aff_mat.*params.mask;
end

if isfield(params, 'labels')
    aff_mat = addLabels2aff(aff_mat, params);
    
end