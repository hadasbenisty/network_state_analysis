function [R2_full, output_vars] = train_full_model(X, params)
    
N = size(X, 1);

switch params.model
    
    case 'lasso'
        for n = 1:N
            predictors = X(setdiff(1:N, n), :);
            
            [B,FitInfo(n)] = lasso(predictors.',X(n, :),'CV',params.folds_num,  'NumLambda',10); %#ok<*AGROW>
            b(:, n) = B(:, FitInfo(n).IndexMinMSE);
            x0(n) = FitInfo(n).Intercept(FitInfo(n).IndexMinMSE);
%             R2_full(n) = var(predictors.'*b(:, n) + x0(n))/var(X(n, :));
SSres = mean((predictors.'*b(:, n) + x0(n) - X(n, :)').^2);
                SStot = var(X(n, :));                
                R2_full(n) = 1 - SSres/SStot;
        end
        output_vars.b = b;
        output_vars.FitInfo = FitInfo;
        output_vars.x0 = x0;
    otherwise
        error('No other options than lasso');
end
end