function visualize_all_experiments(all_results, operating_conditions)
% VISUALIZE_ALL_EXPERIMENTS Creates comprehensive visualization comparing all experiments
%
% Inputs:
%   all_results - Cell array of results structures from each experiment
%   operating_conditions - Matrix [n_exp x 3] with columns: [Pressure(bar), Temp(K), Flow(kg/s)]
%
% Example:
%   visualize_all_experiments(all_results, op_conditions);

    n_exp = length(all_results);

    % Extract data from all experiments
    IF_norms_all = [];
    condition_numbers = zeros(n_exp, 1);
    eigenvalue_ratios = zeros(n_exp, 1);
    max_influences = zeros(n_exp, 1);
    mean_influences = zeros(n_exp, 1);

    % Get number of observations from first experiment
    n_obs = length(all_results{1}.IF_norms);

    % Initialize matrix to store all IF norms
    IF_norms_matrix = zeros(n_obs, n_exp);

    for i = 1:n_exp
        IF_norms_matrix(:, i) = all_results{i}.IF_norms;
        condition_numbers(i) = all_results{i}.condition_number;
        eigenvalue_ratios(i) = all_results{i}.eigenvalues(1) / all_results{i}.eigenvalues(2);
        max_influences(i) = max(all_results{i}.IF_norms);
        mean_influences(i) = mean(all_results{i}.IF_norms);
    end

    % Create comprehensive figure
    figure('Position', [50, 50, 1800, 1200]);

    %% Plot 1: Heatmap of influence norms across all experiments
    subplot(3, 4, 1);
    imagesc(IF_norms_matrix');
    colorbar;
    xlabel('Observation Index');
    ylabel('Experiment Number');
    title('Influence Norms Heatmap (All Experiments)');
    colormap(jet);

    %% Plot 2: Line plot of influence patterns
    subplot(3, 4, 2);
    plot(IF_norms_matrix, 'LineWidth', 1.5);
    xlabel('Observation Index');
    ylabel('||IF_i||');
    title('Influence Patterns Across Experiments');
    legend(arrayfun(@(x) sprintf('Exp %d', x), 1:n_exp, 'UniformOutput', false), ...
        'Location', 'eastoutside', 'FontSize', 7);
    grid on;

    %% Plot 3: Maximum influence per experiment
    subplot(3, 4, 3);
    bar(max_influences);
    xlabel('Experiment Number');
    ylabel('Max ||IF_i||');
    title('Maximum Influence per Experiment');
    grid on;

    %% Plot 4: Mean influence per experiment
    subplot(3, 4, 4);
    bar(mean_influences);
    xlabel('Experiment Number');
    ylabel('Mean ||IF_i||');
    title('Average Influence per Experiment');
    grid on;

    %% Plot 5: Condition numbers across experiments
    subplot(3, 4, 5);
    bar(condition_numbers);
    xlabel('Experiment Number');
    ylabel('Condition Number');
    title('Hessian Condition Numbers');
    grid on;

    %% Plot 6: Eigenvalue ratios
    subplot(3, 4, 6);
    bar(eigenvalue_ratios);
    xlabel('Experiment Number');
    ylabel('\lambda_{max} / \lambda_{min}');
    title('Eigenvalue Ratios');
    grid on;

    %% Plot 7: Influence vs Pressure
    if nargin > 1 && ~isempty(operating_conditions)
        subplot(3, 4, 7);
        scatter(operating_conditions(:, 1), max_influences, 100, 'filled');
        xlabel('Pressure (bar)');
        ylabel('Max ||IF_i||');
        title('Influence vs Pressure');
        grid on;

        % Add experiment labels
        for i = 1:n_exp
            text(operating_conditions(i, 1), max_influences(i), ...
                sprintf(' %d', i), 'FontSize', 8);
        end
    end

    %% Plot 8: Influence vs Temperature
    if nargin > 1 && ~isempty(operating_conditions)
        subplot(3, 4, 8);
        scatter(operating_conditions(:, 2) - 273, max_influences, 100, 'filled');
        xlabel('Temperature (°C)');
        ylabel('Max ||IF_i||');
        title('Influence vs Temperature');
        grid on;

        for i = 1:n_exp
            text(operating_conditions(i, 2) - 273, max_influences(i), ...
                sprintf(' %d', i), 'FontSize', 8);
        end
    end

    %% Plot 9: Influence vs Flow Rate
    if nargin > 1 && ~isempty(operating_conditions)
        subplot(3, 4, 9);
        scatter(operating_conditions(:, 3), max_influences, 100, 'filled');
        xlabel('Flow Rate (kg/s)');
        ylabel('Max ||IF_i||');
        title('Influence vs Flow Rate');
        grid on;

        for i = 1:n_exp
            text(operating_conditions(i, 3), max_influences(i), ...
                sprintf(' %d', i), 'FontSize', 8);
        end
    end

    %% Plot 10: 3D scatter of operating conditions colored by max influence
    if nargin > 1 && ~isempty(operating_conditions)
        subplot(3, 4, 10);
        scatter3(operating_conditions(:, 1), operating_conditions(:, 2) - 273, ...
            operating_conditions(:, 3), 150, max_influences, 'filled');
        xlabel('Pressure (bar)');
        ylabel('Temperature (°C)');
        zlabel('Flow Rate (kg/s)');
        title('Operating Conditions Space');
        colorbar;
        grid on;
        view(45, 30);

        % Add experiment labels
        for i = 1:n_exp
            text(operating_conditions(i, 1), operating_conditions(i, 2) - 273, ...
                operating_conditions(i, 3), sprintf(' %d', i), 'FontSize', 8);
        end
    end

    %% Plot 11: Box plot of influences across experiments
    subplot(3, 4, 11);
    boxplot(IF_norms_matrix, 'Labels', 1:n_exp);
    xlabel('Experiment Number');
    ylabel('||IF_i||');
    title('Distribution of Influences');
    grid on;

    %% Plot 12: Most influential observation per experiment
    subplot(3, 4, 12);
    most_influential_obs = zeros(n_exp, 1);
    for i = 1:n_exp
        most_influential_obs(i) = all_results{i}.most_influential_obs(1);
    end
    bar(most_influential_obs);
    xlabel('Experiment Number');
    ylabel('Observation Index');
    title('Most Influential Observation per Experiment');
    grid on;
    ylim([0, n_obs + 1]);

    sgtitle('Comprehensive Influence Analysis Across All Experiments', ...
        'FontSize', 16, 'FontWeight', 'bold');

    %% Create second figure for parameter-specific analysis
    figure('Position', [100, 100, 1600, 900]);

    % Extract parameter influences (assuming 2 parameters: theta1 and theta2)
    n_params = size(all_results{1}.IF, 1);
    param_names = {'\theta_1 (D_i)', '\theta_2 (\gamma)'};

    for param_idx = 1:n_params
        % Get influences for this parameter across all experiments
        param_IF_matrix = zeros(n_obs, n_exp);
        for i = 1:n_exp
            param_IF_matrix(:, i) = all_results{i}.IF(param_idx, :);
        end

        % Heatmap for this parameter
        subplot(2, 3, (param_idx-1)*3 + 1);
        imagesc(param_IF_matrix');
        colorbar;
        xlabel('Observation Index');
        ylabel('Experiment Number');
        title(sprintf('IF for %s Across Experiments', param_names{param_idx}));
        colormap(jet);

        % Line plot
        subplot(2, 3, (param_idx-1)*3 + 2);
        plot(param_IF_matrix, 'LineWidth', 1.5);
        xlabel('Observation Index');
        ylabel(sprintf('IF(%s)', param_names{param_idx}));
        title(sprintf('Influence on %s', param_names{param_idx}));
        grid on;

        % Summary statistics
        subplot(2, 3, (param_idx-1)*3 + 3);
        mean_param_IF = mean(abs(param_IF_matrix), 1);
        max_param_IF = max(abs(param_IF_matrix), [], 1);

        bar([mean_param_IF; max_param_IF]');
        xlabel('Experiment Number');
        ylabel('Influence Magnitude');
        title(sprintf('Statistics for %s', param_names{param_idx}));
        legend('Mean |IF|', 'Max |IF|', 'Location', 'best');
        grid on;
    end

    sgtitle('Parameter-Specific Influence Analysis', ...
        'FontSize', 16, 'FontWeight', 'bold');

    %% Create third figure for correlation analysis
    if nargin > 1 && ~isempty(operating_conditions)
        figure('Position', [150, 150, 1400, 500]);

        % Calculate correlations
        subplot(1, 3, 1);
        [r_pressure, p_pressure] = corrcoef(operating_conditions(:, 1), max_influences);
        scatter(operating_conditions(:, 1), max_influences, 100, 'filled');
        hold on;
        p = polyfit(operating_conditions(:, 1), max_influences, 1);
        x_fit = linspace(min(operating_conditions(:, 1)), max(operating_conditions(:, 1)), 100);
        plot(x_fit, polyval(p, x_fit), 'r--', 'LineWidth', 2);
        xlabel('Pressure (bar)');
        ylabel('Max ||IF_i||');
        title(sprintf('Pressure vs Influence\n(r = %.3f, p = %.3f)', ...
            r_pressure(1,2), p_pressure(1,2)));
        grid on;

        subplot(1, 3, 2);
        [r_temp, p_temp] = corrcoef(operating_conditions(:, 2), max_influences);
        scatter(operating_conditions(:, 2) - 273, max_influences, 100, 'filled');
        hold on;
        p = polyfit(operating_conditions(:, 2), max_influences, 1);
        x_fit = linspace(min(operating_conditions(:, 2)), max(operating_conditions(:, 2)), 100);
        plot(x_fit - 273, polyval(p, x_fit), 'r--', 'LineWidth', 2);
        xlabel('Temperature (°C)');
        ylabel('Max ||IF_i||');
        title(sprintf('Temperature vs Influence\n(r = %.3f, p = %.3f)', ...
            r_temp(1,2), p_temp(1,2)));
        grid on;

        subplot(1, 3, 3);
        [r_flow, p_flow] = corrcoef(operating_conditions(:, 3), max_influences);
        scatter(operating_conditions(:, 3), max_influences, 100, 'filled');
        hold on;
        p = polyfit(operating_conditions(:, 3), max_influences, 1);
        x_fit = linspace(min(operating_conditions(:, 3)), max(operating_conditions(:, 3)), 100);
        plot(x_fit, polyval(p, x_fit), 'r--', 'LineWidth', 2);
        xlabel('Flow Rate (kg/s)');
        ylabel('Max ||IF_i||');
        title(sprintf('Flow Rate vs Influence\n(r = %.3f, p = %.3f)', ...
            r_flow(1,2), p_flow(1,2)));
        grid on;

        sgtitle('Correlation Between Operating Conditions and Influence', ...
            'FontSize', 14, 'FontWeight', 'bold');
    end

    %% Print summary statistics
    fprintf('\n========================================\n');
    fprintf('  SUMMARY ACROSS ALL EXPERIMENTS\n');
    fprintf('========================================\n\n');

    fprintf('Overall Statistics:\n');
    fprintf('  Mean of max influences: %.4f\n', mean(max_influences));
    fprintf('  Std of max influences:  %.4f\n', std(max_influences));
    fprintf('  Range: [%.4f, %.4f]\n', min(max_influences), max(max_influences));
    fprintf('\n');

    fprintf('Condition Numbers:\n');
    fprintf('  Mean: %.2f\n', mean(condition_numbers));
    fprintf('  Std:  %.2f\n', std(condition_numbers));
    fprintf('  Range: [%.2f, %.2f]\n', min(condition_numbers), max(condition_numbers));
    fprintf('\n');

    if nargin > 1 && ~isempty(operating_conditions)
        fprintf('Correlations with Operating Conditions:\n');
        fprintf('  Pressure:    r = %.3f, p = %.3f\n', r_pressure(1,2), p_pressure(1,2));
        fprintf('  Temperature: r = %.3f, p = %.3f\n', r_temp(1,2), p_temp(1,2));
        fprintf('  Flow Rate:   r = %.3f, p = %.3f\n', r_flow(1,2), p_flow(1,2));
        fprintf('\n');
    end

end
