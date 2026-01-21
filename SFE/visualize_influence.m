function visualize_influence(results, varargin)
% VISUALIZE_INFLUENCE Creates diagnostic plots for influence analysis
%
% Inputs:
%   results - Output structure from compute_influence_diagnostics
%   varargin - Optional name-value pairs:
%       'TopN' - Number of top influential observations to highlight (default: 5)
%       'ParamNames' - Cell array of parameter names (default: {'theta1', 'theta2', ...})

    % Parse optional inputs

    warning('off', 'all')

    p = inputParser;
    addParameter(p, 'TopN', 5, @isnumeric);
    addParameter(p, 'ParamNames', [], @iscell);
    parse(p, varargin{:});
    
    top_n = p.Results.TopN;
    param_names = p.Results.ParamNames;
    
    [n_params, n_obs] = size(results.IF);
    
    % Generate default parameter names if not provided
    if isempty(param_names)
        param_names = cell(n_params, 1);
        for i = 1:n_params
            param_names{i} = sprintf('\\theta_%d', i);
        end
    end
    
    % Create figure with subplots
    figure('Position', [100, 100, 1400, 1000]);
    
    %% Plot 1: Influence function norms
    subplot(3, 3, 1);
    bar(results.IF_norms);
    hold on;
    top_obs = results.most_influential_obs(1:min(top_n, length(results.most_influential_obs)));
    scatter(top_obs, results.IF_norms(top_obs), 100, 'r', 'filled');
    xlabel('Observation');
    ylabel('||IF_i||');
    title('Total Influence Magnitude');
    grid on;
    
    % Highlight top observations
    for i = 1:min(top_n, length(top_obs))
        text(top_obs(i), results.IF_norms(top_obs(i)), ...
            sprintf('  %d', top_obs(i)), 'FontSize', 10);
    end
    
    %% Plot 2: Influence by parameter
    subplot(3, 3, 2);
    plot(results.IF', 'LineWidth', 1.5);
    xlabel('Observation');
    ylabel('IF value');
    title('Influence on Each Parameter');
    legend(param_names, 'Location', 'best');
    grid on;
    
    %% Plot 3: Heatmap of influence
    subplot(3, 3, 3);
    imagesc(results.IF);
    colorbar;
    xlabel('Observation');
    ylabel('Parameter');
    title('Influence Function Heatmap');
    yticks(1:n_params);
    yticklabels(param_names);
    colormap(jet);
    
    %% Plot 4: Parameter sensitivity distribution
    subplot(3, 3, 4);
    histogram(results.parameter_sensitivity, 'BinMethod', 'integers');
    xlabel('Parameter Index');
    ylabel('Frequency');
    title('Which Parameter Most Affected by Each Obs');
    xticks(1:n_params);
    xticklabels(param_names);
    grid on;
    
    %% Plot 5: Eigenvalue spectrum
    subplot(3, 3, 5);
    bar(abs(results.eigenvalues));
    xlabel('Eigenvalue Index');
    ylabel('|Eigenvalue|');
    title(sprintf('Hessian Eigenvalues (Cond: %.2f)', results.condition_number));
    grid on;
    %yscale('log');
    
    %% Plot 6: Variance contributions
    subplot(3, 3, 6);
    total_var_contrib = sum(results.variance_contributions, 2);
    bar(total_var_contrib);
    xlabel('Parameter');
    ylabel('Total Variance Contribution');
    title('Parameter Variance from All Observations');
    xticks(1:n_params);
    xticklabels(param_names);
    grid on;
    
    %% Plot 7: Scatter plot of IF components (for 2D case)
    if n_params == 2
        subplot(3, 3, 7);
        scatter(results.IF(1, :), results.IF(2, :), 50, 1:n_obs, 'filled');
        hold on;
        scatter(results.IF(1, top_obs), results.IF(2, top_obs), ...
            150, 'r', 'LineWidth', 2);
        xlabel(sprintf('IF for %s', param_names{1}));
        ylabel(sprintf('IF for %s', param_names{2}));
        title('2D Influence Space');
        colorbar;
        grid on;
        
        % Add observation labels for top influential
        for i = 1:min(top_n, length(top_obs))
            text(results.IF(1, top_obs(i)), results.IF(2, top_obs(i)), ...
                sprintf('  %d', top_obs(i)), 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
    
    %% Plot 8: Cumulative influence
    subplot(3, 3, 8);
    sorted_norms = sort(results.IF_norms, 'descend');
    cumsum_norms = cumsum(sorted_norms);
    cumsum_normalized = cumsum_norms / cumsum_norms(end);
    plot(1:n_obs, cumsum_normalized, 'LineWidth', 2);
    xlabel('Number of Observations');
    ylabel('Cumulative Proportion of Total Influence');
    title('Cumulative Influence Distribution');
    grid on;
    
    % Add reference line at 80%
    hold on;
    yline(0.8, '--r', '80%');
    n_80 = find(cumsum_normalized >= 0.8, 1);
    xline(n_80, '--r', sprintf('%d obs', n_80));
    
    %% Plot 9: Influence vs observation index
    subplot(3, 3, 9);
    stem(results.IF_norms, 'filled');
    hold on;
    scatter(top_obs, results.IF_norms(top_obs), 100, 'r', 'filled');
    xlabel('Observation');
    ylabel('||IF_i||');
    title('Influence Pattern');
    grid on;
    
    sgtitle('Influence Function Diagnostic Plots', 'FontSize', 14, 'FontWeight', 'bold');
end