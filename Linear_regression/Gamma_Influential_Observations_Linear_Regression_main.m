close all; clc; clear all;

%%
DI = [0.71383013	1.076801997	2.179470155	2.475532632	1.390707877	1.336111172	1.882954204	2.457886055	0.564935512	1.542106938	0.835725102	0.87349666];
GG = [4.229739602	3.091520556	2.359538225	1.132795818	2.204975712	2.739220425	1.868538631	1.69935869	3.452308202	1.995905641	3.012676539	2.596460037];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%%
YY = GG;
name_YY = 'GG';

n = numel(YY);
p = 3;

% Generate design matrix (including intercept)
X = [ones(1,n); RE; FF]';

y = YY';

% Compute diagnostics
diagnostics = regression_influence_diagnostics(X, y);

% Print summary
print_influence_summary(diagnostics);

% Create plots
plot_influence_diagnostics(diagnostics, name_YY);

% Print coefficient estimates
fprintf('=== COEFFICIENT ESTIMATES ===\n');
fprintf('Coefficient estimates:\n');
for j = 1:p
    fprintf('  beta_%d = %.4f\n', j-1, diagnostics.beta_hat(j));
end
fprintf('\n');

% Compare with MATLAB's built-in fitlm (if Statistics Toolbox available)
try
    mdl = fitlm(X(:, 2:end), y);
    fprintf('Comparison with MATLAB fitlm:\n');
    fprintf('  Intercept: %.4f (ours: %.4f)\n', mdl.Coefficients.Estimate(1), diagnostics.beta_hat(1));
    for j = 2:p
        fprintf('  Beta_%d: %.4f (ours: %.4f)\n', j-1, mdl.Coefficients.Estimate(j), diagnostics.beta_hat(j));
    end
catch
    fprintf('Statistics Toolbox not available for comparison.\n');
end

%% Investigate Observation 1
% Look at the raw data for observation 1
fprintf('\n*******************************\n');
fprintf('Observation 1 data:\n');
fprintf('  x1 = %.4f\n', X(1, 2));
fprintf('  x2 = %.4f\n', X(1, 3));
fprintf('  y  = %.4f\n', y(1));
fprintf('  Predicted y = %.4f\n', diagnostics.fitted_values(1));
fprintf('  Residual = %.4f\n', diagnostics.residuals(1));

% Compare to other observations
fprintf('\nComparison to rest of data:\n');
fprintf('  x1 range (others): [%.4f, %.4f]\n', min(X(2:end,2)), max(X(2:end,2)));
fprintf('  x2 range (others): [%.4f, %.4f]\n', min(X(2:end,3)), max(X(2:end,3)));
fprintf('  y range (others): [%.4f, %.4f]\n', min(y(2:end)), max(y(2:end)));

%% Fit Model Without Observation 1
% Remove observation 1
X_reduced = X([2:end], :);
y_reduced = y([2:end]);

% Refit
diagnostics_reduced = regression_influence_diagnostics(X_reduced, y_reduced);
fprintf('\n*******************************\n');
fprintf('\n=== COEFFICIENTS COMPARISON ===\n');
fprintf('               Full Data    Without Obs 1    Change\n');
fprintf('Intercept:     %.4f       %.4f          %.4f\n', ...
        diagnostics.beta_hat(1), diagnostics_reduced.beta_hat(1), ...
        diagnostics.beta_hat(1) - diagnostics_reduced.beta_hat(1));
fprintf('Beta_1:        %.4f       %.4f          %.4f  (%.1f%%)\n', ...
        diagnostics.beta_hat(2), diagnostics_reduced.beta_hat(2), ...
        diagnostics.beta_hat(2) - diagnostics_reduced.beta_hat(2), ...
        100*abs(diagnostics.beta_hat(2) - diagnostics_reduced.beta_hat(2))/abs(diagnostics.beta_hat(2)));
fprintf('Beta_2:        %.4f       %.4f          %.4f  (%.1f%%)\n', ...
        diagnostics.beta_hat(3), diagnostics_reduced.beta_hat(3), ...
        diagnostics.beta_hat(3) - diagnostics_reduced.beta_hat(3), ...
        100*abs(diagnostics.beta_hat(3) - diagnostics_reduced.beta_hat(3))/abs(diagnostics.beta_hat(3)));

%%
% Detailed investigation of DFBETAS-flagged observations
fprintf('\n*******************************\n');
fprintf('═══════════════════════════════════════════════════════\n');
fprintf('DETAILED DFBETAS INVESTIGATION FOR DATASET 2\n');
fprintf('═══════════════════════════════════════════════════════\n\n');

flagged_obs = [4, 8, 9];
threshold_dfbetas = 2/sqrt(diagnostics.n);

for i = 1:length(flagged_obs)
    obs_idx = flagged_obs(i);
    
    fprintf('┌─────────────────────────────────────────────────────┐\n');
    fprintf('│ OBSERVATION %d                                       │\n', obs_idx);
    fprintf('└─────────────────────────────────────────────────────┘\n\n');
    
    % Basic data
    fprintf('Raw Data:\n');
    fprintf('  x₁ = %.4f\n', X(obs_idx, 2));
    fprintf('  x₂ = %.4f\n', X(obs_idx, 3));
    fprintf('  y  = %.4f\n', y(obs_idx));
    
    % Percentiles in the data distribution
    percentile_x1 = 100 * sum(X(:,2) <= X(obs_idx,2)) / diagnostics.n;
    percentile_x2 = 100 * sum(X(:,3) <= X(obs_idx,3)) / diagnostics.n;
    percentile_y = 100 * sum(y <= y(obs_idx)) / diagnostics.n;
    
    fprintf('\nPosition in Data Distribution (percentile):\n');
    fprintf('  x₁: %.1f%% (min=%.4f, median=%.4f, max=%.4f)\n', ...
            percentile_x1, min(X(:,2)), median(X(:,2)), max(X(:,2)));
    fprintf('  x₂: %.1f%% (min=%.4f, median=%.4f, max=%.4f)\n', ...
            percentile_x2, min(X(:,3)), median(X(:,3)), max(X(:,3)));
    fprintf('  y:  %.1f%% (min=%.4f, median=%.4f, max=%.4f)\n', ...
            percentile_y, min(y), median(y), max(y));
    
    % Model fit
    fprintf('\nModel Fit:\n');
    fprintf('  Predicted y = %.4f\n', diagnostics.fitted_values(obs_idx));
    fprintf('  Residual    = %.4f (%.2f SD)\n', ...
            diagnostics.residuals(obs_idx), diagnostics.std_residuals(obs_idx));
    
    % Influence measures
    fprintf('\nInfluence Measures:\n');
    fprintf('  Leverage    = %.4f (threshold: %.4f)\n', ...
            diagnostics.leverage(obs_idx), 2*diagnostics.p/diagnostics.n);
    fprintf('  Cook''s D    = %.4f (threshold: %.4f)\n', ...
            diagnostics.cooks_distance(obs_idx), 4/diagnostics.n);
    fprintf('  DFFITS      = %.4f (threshold: %.4f)\n', ...
            diagnostics.dffits(obs_idx), 2*sqrt(diagnostics.p/diagnostics.n));
    
    % DFBETAS for each coefficient
    fprintf('\nDFBETAS Values:\n');
    coef_names = {'\beta_0 (Intercept)', '\beta_1 (x_1)', '\beta_2 (x_2)'};
    for j = 1:diagnostics.p
        dfbeta_val = diagnostics.dfbetas(obs_idx, j);
        exceeds = abs(dfbeta_val) > threshold_dfbetas;
        marker = '';
        if exceeds
            marker = ' ⚠️ EXCEEDS';
        end
        fprintf('  %-20s: %+.4f (threshold: %.4f)%s\n', ...
                coef_names{j}, dfbeta_val, threshold_dfbetas, marker);
    end
    
    % Distance from centroid
    x_centered = X(obs_idx, 2:end) - mean(X(:, 2:end));
    mahal_dist = sqrt(x_centered * inv(cov(X(:,2:end))) * x_centered');
    fprintf('\nMultivariate Distance:\n');
    fprintf('  Mahalanobis distance from centroid: %.4f\n', mahal_dist);
    
    fprintf('\n');
end

%%
% Compare coefficient estimates with/without each observation
fprintf('\n*******************************\n');
fprintf('═══════════════════════════════════════════════════════\n');
fprintf('COEFFICIENT SENSITIVITY ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════\n\n');

fprintf('Full Model Coefficients:\n');
fprintf('  β₀ = %.4f\n', diagnostics.beta_hat(1));
fprintf('  β₁ = %.4f\n', diagnostics.beta_hat(2));
fprintf('  β₂ = %.4f\n\n', diagnostics.beta_hat(3));

% Create comparison table
fprintf('%-15s %12s %12s %12s %12s %12s %12s\n', ...
        'Model', 'β₀', 'Δβ₀', 'β₁', 'Δβ₁', 'β₂', 'Δβ₂');
fprintf('%-15s %12s %12s %12s %12s %12s %12s\n', ...
        '───────────────', '──────────', '──────────', ...
        '──────────', '──────────', '──────────', '──────────');

for i = 1:length(flagged_obs)
    obs_idx = flagged_obs(i);
    
    % Refit without this observation
    X_reduced = X([1:obs_idx-1, obs_idx+1:end], :);
    y_reduced = y([1:obs_idx-1, obs_idx+1:end]);
    beta_reduced = (X_reduced' * X_reduced) \ (X_reduced' * y_reduced);
    
    % Calculate changes
    delta_beta = diagnostics.beta_hat - beta_reduced;
    
    fprintf('Without Obs %-3d %12.4f %+12.4f %12.4f %+12.4f %12.4f %+12.4f\n', ...
            obs_idx, ...
            beta_reduced(1), delta_beta(1), ...
            beta_reduced(2), delta_beta(2), ...
            beta_reduced(3), delta_beta(3));
end

fprintf('\n');

% Percentage changes
fprintf('Percentage Changes When Removing Each Observation:\n');
fprintf('%-15s %12s %12s %12s\n', 'Remove Obs', 'Δβ₀ (%)', 'Δβ₁ (%)', 'Δβ₂ (%)');
fprintf('%-15s %12s %12s %12s\n', '───────────────', '──────────', '──────────', '──────────');

for i = 1:length(flagged_obs)
    obs_idx = flagged_obs(i);
    
    X_reduced = X([1:obs_idx-1, obs_idx+1:end], :);
    y_reduced = y([1:obs_idx-1, obs_idx+1:end]);
    beta_reduced = (X_reduced' * X_reduced) \ (X_reduced' * y_reduced);
    
    pct_change = 100 * abs(diagnostics.beta_hat - beta_reduced) ./ abs(diagnostics.beta_hat);
    
    fprintf('%-15d %11.2f%% %11.2f%% %11.2f%%\n', ...
            obs_idx, pct_change(1), pct_change(2), pct_change(3));
end

fprintf('\n');

%%