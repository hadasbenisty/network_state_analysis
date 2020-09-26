function [R2_partial, output_vars] = train_partial_model(X, params)

N = size(X, 1);

switch params.model
    
    case 'lasso'
        b = zeros(N-2, N, N);
        for n = 1:N
            for k = 1:N
                if k == n
                    R2_partial(n, k) = 0;
                    continue;
                end
                predictors = X(setdiff(1:N, [n, k]), :);
                
                [B,FitInfo(n, k)] = lasso(predictors.',X(n, :),'CV',params.folds_num,  'NumLambda',10); %#ok<*AGROW>
                b(:, n, k) = B(:, FitInfo(n, k).IndexMinMSE);
                x0(n, k) = FitInfo(n, k).Intercept(FitInfo(n, k).IndexMinMSE);
                SSres = mean((predictors.'*b(:, n, k) + x0(n, k) - X(n, :)').^2);
                SStot = var(X(n, :));                
                R2_partial(n, k) = 1 - SSres/SStot;
            end
        end
        
        output_vars.b = b;
        output_vars.FitInfo = FitInfo;
        output_vars.x0 = x0;
    otherwise
        error('No other options than lasso');
end
end