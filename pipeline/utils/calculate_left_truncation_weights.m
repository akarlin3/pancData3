function weights = calculate_left_truncation_weights(t_start, t_stop)
    % Calculate weights to adjust for left-truncation bias
    % Higher weights for observations with later entry times
    
    % Base weight: inverse probability of surviving to entry
    max_start = max(t_start);
    if max_start > 0
        % Simple approximation: exponential survival to entry
        entry_survival_prob = exp(-0.1 * t_start / max_start);
        weights = 1 ./ (entry_survival_prob + 0.1);  % Avoid extreme weights
    else
        weights = ones(size(t_start));
    end
    
    % Normalize weights
    weights = weights / mean(weights);
    
    % Cap extreme weights
    weights = min(weights, 5);  % Maximum weight = 5
    weights = max(weights, 0.2); % Minimum weight = 0.2
end
