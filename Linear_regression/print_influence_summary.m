function print_influence_summary(diagnostics, threshold_cooks, threshold_dffits)
% PRINT_INFLUENCE_SUMMARY Print summary of influential observations
%
% print_influence_summary(diagnostics)
% print_influence_summary(diagnostics, threshold_cooks, threshold_dffits)
%
% Inputs:
%   diagnostics      - structure from regression_influence_diagnostics
%   threshold_cooks  - (optional) threshold for Cook's distance (default: 4/n)
%   threshold_dffits - (optional) threshold for DFFITS (default: 2*sqrt(p/n))

    n = diagnostics.n;
    p = diagnostics.p;
    
    if nargin < 2
        threshold_cooks = 4/n;
    end
    if nargin < 3
        threshold_dffits = 2 * sqrt(p/n);
    end
    
    fprintf('\n=== INFLUENCE DIAGNOSTICS SUMMARY ===\n');
    fprintf('Sample size (n): %d\n', n);
    fprintf('Number of parameters (p): %d\n', p);
    fprintf('Residual standard error: %.4f\n\n', diagnostics.sigma_hat);
    
    % Cook's Distance
    influential_cooks = find(diagnostics.cooks_distance > threshold_cooks);
    fprintf("Cook's Distance (threshold = %.4f):\n", threshold_cooks);
    if isempty(influential_cooks)
        fprintf('  No influential observations detected.\n');
    else
        fprintf('  %d influential observations:\n', length(influential_cooks));
        for i = 1:length(influential_cooks)
            idx = influential_cooks(i);
            fprintf('    Obs %d: D = %.4f\n', idx, diagnostics.cooks_distance(idx));
        end
    end
    fprintf('\n');
    
    % DFFITS
    influential_dffits = find(abs(diagnostics.dffits) > threshold_dffits);
    fprintf('DFFITS (threshold = %.4f):\n', threshold_dffits);
    if isempty(influential_dffits)
        fprintf('  No influential observations detected.\n');
    else
        fprintf('  %d influential observations:\n', length(influential_dffits));
        for i = 1:length(influential_dffits)
            idx = influential_dffits(i);
            fprintf('    Obs %d: DFFITS = %.4f\n', idx, diagnostics.dffits(idx));
        end
    end
    fprintf('\n');
    
    % High Leverage
    high_leverage = find(diagnostics.leverage > 2*p/n);
    fprintf('High Leverage (threshold = 2p/n = %.4f):\n', 2*p/n);
    if isempty(high_leverage)
        fprintf('  No high leverage observations detected.\n');
    else
        fprintf('  %d high leverage observations:\n', length(high_leverage));
        for i = 1:length(high_leverage)
            idx = high_leverage(i);
            fprintf('    Obs %d: h_ii = %.4f\n', idx, diagnostics.leverage(idx));
        end
    end
    fprintf('\n');
    
    % DFBETAS
    threshold_dfbetas = 2 / sqrt(n);
    fprintf('DFBETAS (threshold = 2/sqrt(n) = %.4f):\n', threshold_dfbetas);
    influential_any_coef = false;
    for j = 1:p
        influential_dfbetas = find(abs(diagnostics.dfbetas(:, j)) > threshold_dfbetas);
        if ~isempty(influential_dfbetas)
            influential_any_coef = true;
            fprintf('  Coefficient %d: %d influential observations\n', j, length(influential_dfbetas));
            for i = 1:min(5, length(influential_dfbetas))  % Show first 5
                idx = influential_dfbetas(i);
                fprintf('    Obs %d: DFBETAS = %.4f\n', idx, diagnostics.dfbetas(idx, j));
            end
            if length(influential_dfbetas) > 5
                fprintf('    ... and %d more\n', length(influential_dfbetas) - 5);
            end
        end
    end
    if ~influential_any_coef
        fprintf('  No influential observations detected for any coefficient.\n');
    end
    fprintf('\n');
end