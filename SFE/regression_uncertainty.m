function [beta_mean, beta_std] = regression_uncertainty(X, y, coef_names)
%REGRESSION_UNCERTAINTY Computes regression coefficients and their uncertainties.
%   [beta_mean, beta_std] = regression_uncertainty(X, y, coef_names)
%   Also validates and plots Gaussian distributions of the coefficient estimates.
%
%   Inputs:
%       X          - Matrix of predictors (n x p), without intercept column
%       y          - Response vector (n x 1)
%       coef_names - Cell array of coefficient names (1 x (p+1)). Optional.
%                    First name should be for the intercept.
%
%   Outputs:
%       beta_mean - Estimated regression coefficients ((p+1) x 1)
%       beta_std  - Standard deviation (standard error) of each coefficient ((p+1) x 1)

    % Add intercept term
    X_aug = [ones(size(X, 1), 1), X];

    % Number of observations and parameters
    [n, p] = size(X_aug);

    % Estimate coefficients (mean)
    beta_mean = X_aug \ y;

    % Residuals and variance of the errors
    residuals = y - X_aug * beta_mean;
    sigma2 = sum(residuals .^ 2) / (n - p);

    % Covariance matrix of coefficients
    XtX = X_aug' * X_aug;
    cov_beta = sigma2 * inv(XtX);

    % === Variance-Covariance Matrix Validation ===
    % Check symmetry
    if ~isequal(cov_beta, cov_beta')
        warning('Covariance matrix is not symmetric.');
    end

    % Check positive semi-definiteness
    eig_vals = eig(cov_beta);
    if any(eig_vals < -1e-10)  % allow small negative due to numerical precision
        warning('Covariance matrix is not positive semi-definite.');
    end

    % Standard deviation (standard error) of coefficients
    beta_std = sqrt(diag(cov_beta));

    cov_beta

    % Default coefficient names if none provided
    if nargin < 3 || isempty(coef_names)
        coef_names = ["Intercept"];
        for i = 1:p-1
            coef_names(end+1) = "X" + i;
        end
    end

    % Plot Gaussian for each coefficient
    figure;
    num_params = length(beta_mean);
    for i = 1:num_params
        subplot(ceil(num_params/2), 2, i);
        mu = beta_mean(i);
        sigma = beta_std(i);
        x_range = linspace(mu - 4*sigma, mu + 4*sigma, 200);
        y_gauss = normpdf(x_range, mu, sigma);
        plot(x_range, y_gauss, 'LineWidth', 2);
        if iscell(coef_names)
            name = coef_names{i};
        else
            name = coef_names(i);
        end
        title(sprintf('%s (%.2f ± %.2f)', name, mu, sigma));
        xlabel('\beta');
        ylabel('PDF');
        grid on;
    end
    sgtitle('Gaussian Distribution of Regression Coefficients');
end
