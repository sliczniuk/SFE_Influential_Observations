function diagnostics = regression_influence_diagnostics(X, y)
% REGRESSION_INFLUENCE_DIAGNOSTICS Compute influence diagnostics for linear regression
%
% diagnostics = regression_influence_diagnostics(X, y)
%
% Inputs:
%   X - n x p design matrix (including intercept column if desired)
%   y - n x 1 response vector
%
% Outputs:
%   diagnostics - structure containing:
%       .beta_hat       - p x 1 coefficient estimates
%       .residuals      - n x 1 residuals
%       .fitted_values  - n x 1 fitted values
%       .sigma_hat      - residual standard error
%       .leverage       - n x 1 leverage values (diagonal of hat matrix)
%       .cooks_distance - n x 1 Cook's distances
%       .dffits         - n x 1 DFFITS values
%       .dfbetas        - n x p DFBETAS matrix
%       .std_residuals  - n x 1 standardized residuals
%       .n              - sample size
%       .p              - number of parameters

    % Get dimensions
    [n, p] = size(X);
    
    % Compute OLS estimates
    XtX = X' * X;
    XtX_inv = inv(XtX);
    beta_hat = XtX_inv * (X' * y);
    
    % Fitted values and residuals
    fitted_values = X * beta_hat;
    residuals = y - fitted_values;
    
    % Residual sum of squares and standard error
    RSS = sum(residuals.^2);
    sigma_hat = sqrt(RSS / (n - p));
    
    % Leverage (diagonal of hat matrix)
    % H = X * (X'X)^{-1} * X'
    % h_ii = X_i * (X'X)^{-1} * X_i'
    leverage = zeros(n, 1);
    for i = 1:n
        leverage(i) = X(i, :) * XtX_inv * X(i, :)';
    end
    
    % Standardized residuals
    % r_i = e_i / (sigma * sqrt(1 - h_ii))
    std_residuals = residuals ./ (sigma_hat * sqrt(1 - leverage));
    
    % Cook's Distance
    % D_i = (r_i^2 / p) * (h_ii / (1 - h_ii))
    cooks_distance = (std_residuals.^2 / p) .* (leverage ./ (1 - leverage));
    
    % DFFITS
    % DFFITS_i = r_i * sqrt(h_ii / (1 - h_ii))
    dffits = std_residuals .* sqrt(leverage ./ (1 - leverage));
    
    % DFBETAS
    % DFBETAS_ij = (beta_j - beta_j(i)) / (sigma_(i) * sqrt((X'X)^{-1}_{jj}))
    % We can compute this efficiently using the Sherman-Morrison result
    dfbetas = zeros(n, p);
    
    for i = 1:n
        % Change in coefficients when removing observation i
        % beta_hat - beta_hat_(i) = (X'X)^{-1} * X_i' * e_i / (1 - h_ii)
        delta_beta = XtX_inv * X(i, :)' * residuals(i) / (1 - leverage(i));
        
        % Leave-one-out residual standard error
        % This is a more accurate computation
        sigma_hat_i = sqrt((RSS - residuals(i)^2 / (1 - leverage(i))) / (n - p - 1));
        
        % Standardize by coefficient standard errors
        for j = 1:p
            dfbetas(i, j) = delta_beta(j) / (sigma_hat_i * sqrt(XtX_inv(j, j)));
        end
    end
    
    % Store results in structure
    diagnostics.beta_hat = beta_hat;
    diagnostics.residuals = residuals;
    diagnostics.fitted_values = fitted_values;
    diagnostics.sigma_hat = sigma_hat;
    diagnostics.leverage = leverage;
    diagnostics.cooks_distance = cooks_distance;
    diagnostics.dffits = dffits;
    diagnostics.dfbetas = dfbetas;
    diagnostics.std_residuals = std_residuals;
    diagnostics.n = n;
    diagnostics.p = p;
end