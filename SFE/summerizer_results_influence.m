function [results] = summerizer_results_influence(H, jacobians)
     % Print input matrices
    fprintf('\n========================================\n');
    fprintf('         INPUT MATRICES\n');
    fprintf('========================================\n\n');
    
    % Print Hessian
    fprintf('HESSIAN MATRIX (H):\n');
    fprintf('------------------\n');
    [n_params, ~] = size(H);
       
    % Print with proper formatting
    for i = 1:n_params
        fprintf('  ');
        for j = 1:n_params
            fprintf('%12.4f', H(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
    
    % Compute influence diagnostics
    fprintf('========================================\n');
    fprintf('    INFLUENCE FUNCTION ANALYSIS\n');
    fprintf('========================================\n\n');
    
    results = compute_influence_diagnostics(H, jacobians);
    
    % Display results
    fprintf('HESSIAN DIAGNOSTICS:\n');
    fprintf('-------------------\n');
    fprintf('Condition Number: %.2f\n', results.condition_number);
    fprintf('Eigenvalues: [%.2f, %.2f]\n', results.eigenvalues(1), results.eigenvalues(2));
    fprintf('Eigenvalue ratio: %.2f (larger/smaller)\n', ...
        results.eigenvalues(1)/results.eigenvalues(2));
    fprintf('\n');
    
    fprintf('Eigenvectors:\n');
    for i = 1:n_params
        fprintf('  v_%d = [', i);
        for j = 1:n_params
            fprintf('%7.4f', results.eigenvectors(j,i));
            if j < n_params
                fprintf(', ');
            end
        end
        fprintf(']\n');
    end
    fprintf('\n');
    
    % Display full table
    fprintf('========================================\n');
    fprintf('      COMPLETE DATA TABLE\n');
    fprintf('========================================\n\n');
    disp(results.summary_table);
end