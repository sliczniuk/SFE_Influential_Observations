function plot_influence_diagnostics(diagnostics, name_YY, observation_labels)
% PLOT_INFLUENCE_DIAGNOSTICS Create diagnostic plots for influence measures
%
% plot_influence_diagnostics(diagnostics)
% plot_influence_diagnostics(diagnostics, observation_labels)
%
% Inputs:
%   diagnostics        - structure from regression_influence_diagnostics
%   name_YY            - name used to save files
%   observation_labels - (optional) cell array of labels for observations

    if nargin < 3
        observation_labels = cellstr(num2str((1:diagnostics.n)'));
    end
    
    n = diagnostics.n;
    p = diagnostics.p;
    
    % Create figure with subplots
    %figure('Position', [100, 100, 1200, 900]);
    
    % --- Plot 1: Cook's Distance ---
    %subplot(2, 3, 1);
    figure(); set(gca,'FontSize',12);
    stem(1:n, diagnostics.cooks_distance, 'filled');
    hold on;
    % Add threshold lines
    yline(1, 'r--', 'LineWidth', 1.5, 'Label', 'D = 1');
    yline(4/n, 'g--', 'LineWidth', 1.5, 'Label', sprintf('D = 4/n = %.3f', 4/n));
    hold off;
    xlabel('Observation Index');
    ylabel("Cook's Distance");
    title("Cook's Distance");
    grid on;
    
    % Label influential points
    influential_cooks = find(diagnostics.cooks_distance > 4/n);
    if ~isempty(influential_cooks)
        hold on;
        text(influential_cooks, diagnostics.cooks_distance(influential_cooks), ...
             observation_labels(influential_cooks), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_Cooks_D.png'], "Resolution",500); close all;

    % --- Plot 2: DFFITS ---
    %subplot(2, 3, 2);
    figure(); set(gca,'FontSize',12);
    stem(1:n, diagnostics.dffits, 'filled');
    hold on;
    % Add threshold lines
    threshold_dffits = 2 * sqrt(p/n);
    yline(threshold_dffits, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('2*sqrt(p/n) = %.3f', threshold_dffits));
    yline(-threshold_dffits, 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Observation Index');
    ylabel('DFFITS');
    title('DFFITS');
    grid on;
    
    % Label influential points
    influential_dffits = find(abs(diagnostics.dffits) > threshold_dffits);
    if ~isempty(influential_dffits)
        hold on;
        text(influential_dffits, diagnostics.dffits(influential_dffits), ...
             observation_labels(influential_dffits), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_DFFITS.png'], "Resolution",500); close all;

    % --- Plot 3: Leverage ---
    %subplot(2, 3, 3);
    figure(); set(gca,'FontSize',12);
    stem(1:n, diagnostics.leverage, 'filled');
    hold on;
    % Add threshold lines
    yline(2*p/n, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('2p/n = %.3f', 2*p/n));
    yline(3*p/n, 'g--', 'LineWidth', 1.5, ...
          'Label', sprintf('3p/n = %.3f', 3*p/n));
    hold off;
    xlabel('Observation Index');
    ylabel('Leverage ($h_{ii}$)');
    title('Leverage Values');
    grid on;
    
    % Label high leverage points
    high_leverage = find(diagnostics.leverage > 2*p/n);
    if ~isempty(high_leverage)
        hold on;
        text(high_leverage, diagnostics.leverage(high_leverage), ...
             observation_labels(high_leverage), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_Leverage.png'], "Resolution",500); close all;

    % --- Plot 4: Residuals vs Fitted ---
    %subplot(2, 3, 4);
    figure(); set(gca,'FontSize',12);
    scatter(diagnostics.fitted_values, diagnostics.std_residuals, 'filled');
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    yline(2, 'r--', 'LineWidth', 1);
    yline(-2, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('Fitted Values');
    ylabel('Standardized Residuals');
    title('Residuals vs Fitted Values');
    grid on;
    
    % Label outliers
    outliers = find(abs(diagnostics.std_residuals) > 2);
    if ~isempty(outliers)
        hold on;
        text(diagnostics.fitted_values(outliers), diagnostics.std_residuals(outliers), ...
             observation_labels(outliers), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_Residuals.png'], "Resolution",500); close all;

    % --- Plot 5: Cook's Distance vs Leverage ---
    %subplot(2, 3, 5);
    figure(); set(gca,'FontSize',12);
    scatter(diagnostics.leverage, diagnostics.cooks_distance, 'filled');
    xlabel('Leverage ($h_{ii}$)');
    ylabel("Cook's Distance");
    title("Cook's Distance vs Leverage");
    grid on;
    
    % Add contour lines for constant standardized residuals
    hold on;
    h_range = linspace(0, max(diagnostics.leverage), 100);
    for r_val = [1, 2, 3]
        D_contour = (r_val^2 / p) * (h_range ./ (1 - h_range));
        plot(h_range, D_contour, '--', 'LineWidth', 1, ...
             'DisplayName', sprintf('|r| = %d', r_val));
    end
    legend('Location', 'northwest');
    hold off;
    
    % Label influential points
    if ~isempty(influential_cooks)
        hold on;
        text(diagnostics.leverage(influential_cooks), ...
             diagnostics.cooks_distance(influential_cooks), ...
             observation_labels(influential_cooks), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_CD_VS_L.png'], "Resolution",500); close all;
    
    % --- Plot 6: DFBETAS for first coefficient (often intercept) ---
    %subplot(2, 3, 6);
    figure(); set(gca,'FontSize',12);
    stem(1:n, diagnostics.dfbetas(:, 1), 'filled');
    hold on;
    threshold_dfbetas = 2 / sqrt(n);
    yline(threshold_dfbetas, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('2/sqrt(n) = %.3f', threshold_dfbetas));
    yline(-threshold_dfbetas, 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Observation Index');
    ylabel('DFBETAS');
    title('DFBETAS for Coefficient 1');
    grid on;
    
    % Label influential points
    influential_dfbetas = find(abs(diagnostics.dfbetas(:, 1)) > threshold_dfbetas);
    if ~isempty(influential_dfbetas)
        hold on;
        text(influential_dfbetas, diagnostics.dfbetas(influential_dfbetas, 1), ...
             observation_labels(influential_dfbetas), 'VerticalAlignment', 'bottom');
        hold off;
    end
    exportgraphics(figure(1),[name_YY,'_DFBETAS.png'], "Resolution",500); close all;
end