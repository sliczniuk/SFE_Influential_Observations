close all; clc; clear all;

%%
DI = [0.71383013	1.076801997	2.179470155	2.475532632	1.390707877	1.336111172	1.882954204	2.457886055	0.564935512	1.542106938	0.835725102	0.87349666];
GG = [4.229739602	3.091520556	2.359538225	1.132795818	2.204975712	2.739220425	1.868538631	1.69935869	3.452308202	1.995905641	3.012676539	2.596460037];
RE = [0.4632, 0.3783, 0.3029, 0.2619, 0.3579, 0.3140, 0.2635, 0.2323, 0.1787, 0.1160, 0.1889, 0.1512];
FF = [6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 6.67, 3.33, 3.33, 3.33, 3.33];

%%
YY = DI;
name_YY = 'DI';

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
% Create a comparison plot
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
bar([diagnostics.beta_hat, diagnostics_reduced.beta_hat]);
legend('With Obs 1', 'Without Obs 1');
ylabel('Coefficient Value');
xlabel('Parameter');
set(gca, 'XTickLabel', {'\beta_0', '\beta_1', '\beta_2'});
title('Coefficient Estimates');
grid on;

subplot(1, 3, 2);
scatter(1:12, diagnostics.cooks_distance, 100, 'filled');
hold on;
scatter(1, diagnostics.cooks_distance(1), 150, 'r', 'filled');
yline(4/12, 'r--', 'Threshold');
hold off;
xlabel('Observation');
ylabel('Cook''s Distance');
title('Influence of Each Observation');
grid on;

subplot(1, 3, 3);
scatter(diagnostics.leverage, abs(diagnostics.std_residuals), 100, 'filled');
hold on;
scatter(diagnostics.leverage(1), abs(diagnostics.std_residuals(1)), 150, 'r', 'filled');
text(diagnostics.leverage(1), abs(diagnostics.std_residuals(1)), ' Obs 1', 'FontSize', 12);
hold off;
xlabel('Leverage');
ylabel('|Standardized Residual|');
title('Leverage vs Residual');
grid on;

%%
% Check Observation 10
fprintf('\n*******************************\n');
fprintf('\nObservation 10:\n');
fprintf('  Cook''s D: %.4f\n', diagnostics.cooks_distance(10));
fprintf('  DFFITS: %.4f\n', diagnostics.dffits(10));
fprintf('  Leverage: %.4f\n', diagnostics.leverage(10));


%% Bootstrap Confidence Intervals 
% Bootstrap to assess uncertainty
n_boot = 1000;
beta1_boot = zeros(n_boot, 1);

rng(42);
for i = 1:n_boot
    % Resample with replacement
    idx = randi(12, 12, 1);
    X_boot = X(idx, :);
    y_boot = y(idx);
    
    % Fit model
    beta_boot = (X_boot' * X_boot) \ (X_boot' * y_boot);
    beta1_boot(i) = beta_boot(2);
end

fprintf('\n=== BOOTSTRAP RESULTS (β₁) ===\n');
fprintf('OLS estimate:      %.4f\n', diagnostics.beta_hat(2));
fprintf('Bootstrap mean:    %.4f\n', mean(beta1_boot));
fprintf('Bootstrap median:  %.4f\n', median(beta1_boot));
fprintf('Bootstrap 95%% CI:  [%.4f, %.4f]\n', prctile(beta1_boot, [2.5, 97.5]));
fprintf('Bootstrap SD:      %.4f\n', std(beta1_boot));