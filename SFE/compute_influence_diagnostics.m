function results = compute_influence_diagnostics(H, jacobians)
% COMPUTE_INFLUENCE_DIAGNOSTICS Computes influence functions and diagnostics
%
% Inputs:
%   H          - Hessian matrix (p x p) at the MLE/estimate
%   jacobians  - Matrix of gradients (p x n) where each column is the 
%                gradient contribution from observation i
%
% Outputs:
%   results    - Structure containing:
%       .IF              - Influence functions (p x n matrix)
%       .IF_norms        - Total influence magnitude per observation
%       .parameter_sensitivity - Which parameter most affected per obs
%       .max_influence_values  - Maximum influence value per obs
%       .most_influential_obs  - Indices of top influential observations
%       .variance_contributions - Contribution to parameter variance
%       .H_inv           - Inverse Hessian
%       .eigenvalues     - Eigenvalues of H
%       .eigenvectors    - Eigenvectors of H
%       .condition_number - Condition number of H
%       .summary_table   - Formatted table of results

    % Get dimensions
    [p, n] = size(jacobians);
    
    % Validate inputs
    if size(H, 1) ~= p || size(H, 2) ~= p
        error('Hessian dimensions must match number of parameters');
    end
    
    %% 1. Compute inverse Hessian
    H_inv = inv(H);
    results.H_inv = H_inv;
    
    %% 2. Eigenvalue decomposition of Hessian
    [V, D] = eig(H);
    eigenvalues = diag(D);
    [eigenvalues_sorted, idx] = sort(abs(eigenvalues), 'descend');
    eigenvectors_sorted = V(:, idx);
    
    results.eigenvalues = eigenvalues_sorted;
    results.eigenvectors = eigenvectors_sorted;
    results.condition_number = eigenvalues_sorted(1) / eigenvalues_sorted(end);
    
    %% 3. Compute influence functions
    % IF_i = -H^{-1} * grad_i for each observation i
    IF = -H_inv * jacobians;
    results.IF = IF;
    
    %% 4. Compute influence norms (total magnitude)
    IF_norms = sqrt(sum(IF.^2, 1));
    results.IF_norms = IF_norms';
    
    %% 5. Identify which parameter is most affected by each observation
    [max_influence_values, parameter_sensitivity] = max(abs(IF), [], 1);
    results.parameter_sensitivity = parameter_sensitivity';
    results.max_influence_values = max_influence_values';
    
    %% 6. Identify most influential observations
    [~, sorted_idx] = sort(IF_norms, 'descend');
    results.most_influential_obs = sorted_idx';
    
    %% 7. Variance contribution analysis
    % Contribution of each observation to parameter variance
    variance_contributions = IF.^2;
    % Normalize to get proportions
    total_var = sum(variance_contributions, 2);
    variance_props = variance_contributions ./ total_var;
    results.variance_contributions = variance_contributions;
    results.variance_proportions = variance_props;
    
    %% 8. Parameter-level statistics
    results.parameter_stats = struct();
    for j = 1:p
        results.parameter_stats(j).param_index = j;
        results.parameter_stats(j).total_influence = norm(IF(j, :));
        results.parameter_stats(j).mean_abs_influence = mean(abs(IF(j, :)));
        results.parameter_stats(j).max_abs_influence = max(abs(IF(j, :)));
        results.parameter_stats(j).most_influential_obs = ...
            find(abs(IF(j, :)) == results.parameter_stats(j).max_abs_influence, 1);
    end
    
    %% 9. Create summary table
    summary_table = table();
    summary_table.Observation = (1:n)';
    
    % Add gradient components
    for j = 1:p
        summary_table.(sprintf('Grad_theta%d', j)) = jacobians(j, :)';
    end
    
    % Add influence function components
    for j = 1:p
        summary_table.(sprintf('IF_theta%d', j)) = IF(j, :)';
    end
    
    summary_table.IF_Norm = IF_norms';
    summary_table.Most_Affected_Param = parameter_sensitivity';
    summary_table.Max_Influence_Value = max_influence_values';
    
    results.summary_table = summary_table;
    
    %% 10. Additional diagnostics
    results.diagnostics = struct();
    
    % Effective sample size (based on influence)
    % Obs with high influence effectively count as multiple observations
    influence_weights = 1 + IF_norms.^2;
    results.diagnostics.effective_n = sum(influence_weights);
    
    % Maximum leverage (largest eigenvalue of each observation's influence)
    results.diagnostics.max_leverage_obs = sorted_idx(1);
    results.diagnostics.max_leverage_value = IF_norms(sorted_idx(1));
    
    % Directional influence (which eigenvector direction most influenced)
    influence_in_eigendirections = eigenvectors_sorted' * IF;
    [~, dominant_direction_per_obs] = max(abs(influence_in_eigendirections), [], 1);
    results.diagnostics.dominant_eigendirection = dominant_direction_per_obs';
    
end