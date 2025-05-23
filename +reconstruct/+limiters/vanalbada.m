function limited_slope = vanalbada(delta_minus, delta_plus, cfg)

    % VANALBADA Van Albada slope limiter for TVD schemes
    %
    % The Van Albada limiter is designed to be differentiable and symmetric,
    % providing smooth transitions.
    % Author: Dr. Denys Dutykh
    %
    % Syntax:
    %   limited_slope = vanalbada(delta_minus, delta_plus)
    % 
    % Inputs:
    %   delta_minus - Backward difference (q_i - q_{i-1})
    %   delta_plus  - Forward difference (q_{i+1} - q_i)
    % 
    % Output:
    %   limited_slope - Limited slope

    % Use global epsilon from config for division-by-zero protection
    epsilon = cfg.numerics.epsilon;

    % Check if slopes have different signs
    if delta_minus * delta_plus <= 0
        % Opposite signs or one is zero - return zero
        limited_slope = 0;
    else
        % Same sign - use van Albada formula
        numerator = (delta_plus^2 * delta_minus + delta_minus^2 * delta_plus);
        denominator = (delta_plus^2 + delta_minus^2 + epsilon);
        limited_slope = numerator / denominator;
    end

end