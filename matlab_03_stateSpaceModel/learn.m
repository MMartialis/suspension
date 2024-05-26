function [best_mutations, best_outs, bests_error] = learn(time, u, params, out, var, G, P, p, errors, delta_t)

valve = u(:, 1); % vector of valve signal for each time sample
in = [u(:, 2), u(:, 3)]; % matrix of input position and velocity for each time sample

params_b = params; % Initialize the parameters
errors_b = errors; % Initialize the error
outs_b = simulateModel(params_b); % Initialize the output
threshold = mean(repmat(out, 1, P)); % Threshold for the output


for gen = 1:G
    tic
    params_mut = normrnd(params_b, var*(abs(params_b))); % Generate random values for the parameters
    outs_g = simulateModel(params_mut); % Simulate the model
    toc
    disp(['simulation ' gen ' complete'])
    % find a scaling
    errors = sum((outs_g - out).^2, 1); % Compute the error
    improvement = (errors < errors_b) & (max(abs(outs_g)) > 0); % Find the improvements
    if any(improvement) % Update the parameters and error
        params_b(improvement, :) = params_mut(improvement, :);
        errors_b(improvement) = errors(improvement);
        outs_b(:, improvement) = outs_g(:, improvement);
    end
end

% Find the best mutations
[sorted_errors, sorted_indices] = sort(errors_b); % Sort the errors in ascending order
best_mutations = params_b(sorted_indices(1:p), :);
best_outs = outs_b(:, sorted_indices(1:p));
bests_error = sorted_errors(1:p);

    function y = simulateModel(params)
        m_coeffs = params(:, 2:2);
        b_coeffs = params(:, 4:5);
        k_coeffs = params(:, 7:7);

        poly_vals = zeros(numel(time), 3*P);
        
        poly_vals(:, 1:3:end) = polyvalParallel(m_coeffs, valve); % matrix of masses for each specimen and time sample
        poly_vals(:, 2:3:end) = polyvalParallel(b_coeffs, valve); % matrix of damping coefficients for each specimen and time sample
        poly_vals(:, 3:3:end) = polyvalParallel(k_coeffs, valve); % matrix of spring constants for each specimen and time sample
        % m = polyvalParallel(m_coeffs, valve); % matrix of masses for each specimen and time sample
        % b = polyvalParallel(b_coeffs, valve); % matrix of damping coefficients for each specimen and time sample
        % k = polyvalParallel(k_coeffs, valve); % matrix of spring constants for each specimen and time sample
        
        % options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); % Set the options for the ODE solvertidal
        
        options = odeset('RelTol',1e-2);
        [~, output] = ode45(@(t, x) stateEquationsTimeVarying(t, x, poly_vals), time, repmat([0; out(1)],1, P), options); % Solve the ODE
        y = output(:, 2:2:end)*delta_t; % Compute the output

    end

    function dxdt = stateEquationsTimeVarying(t, x, poly_vals)
        % Interpolate the parameter values at time t
        % p1 = interp1(time, p1_vals, t, 'linear');
        % p2 = interp1(time, p2_vals, t, 'linear');
        % p3 = interp1(time, p3_vals, t, 'linear');
        % u = interp1(time, in(:, 1), t, 'linear');
        % du = interp1(time, in(:, 2), t, 'linear');

        % p1 = p1_vals(max(1,round(t/delta_t)), :);
        % p2 = p2_vals(max(1,round(t/delta_t)), :);
        % p3 = p3_vals(max(1,round(t/delta_t)), :);
        % u = in(max(1, round(t/delta_t)), 1);
        % du = in(max(1, round(t/delta_t)), 2);
        
        % p_t = poly_vals(max(1,round(t/delta_t)), :);
        % u_t = in(max(1, round(t/delta_t)), :);
        % [p1, p2, p3] = deal(p_t(1), p_t(2), p_t(3));
        % [u, du] = deal(u_t(1), u_t(2));

        index = max(1, round(t/delta_t));
        [p1, p2, p3, u, du] = deal(poly_vals(index, 1), poly_vals(index, 2), poly_vals(index, 3), in(index, 1), in(index, 2));


        % Calculate the derivatives
        dxdt = zeros(numel(x),1);
        dxdt(1:2:end) = x(2:2:end);
        dxdt(2:2:end) = (p2 / p1) * du + (p3 / p1) * u - (p2 / p1) * x(2) - (p3 / p1) * x(1);
    end

end

function y = polyvalParallel(coeffs, x)
% coeffs is an NxM matrix where N is the number of sets of coefficients
% and M is the polynomial order + 1.
% x is a vector of input values.
numVals = size(x, 1);
polyOrder = size(coeffs, 2) - 1;

% Preallocate matrix for powers of x
xPowers = ones(numVals, polyOrder + 1);

% Compute powers iteratively
if polyOrder > 1
    xPowers(:, end-1) = x;
    for i = polyOrder-1:1
        xPowers(:, i) = xPowers(:, i+1) .* x;
    end
else
    xPowers(:, 1) = x;
end
% Add the column for x^0 (which is 1)
xPowers(:, end) = 1;
% Perform matrix multiplication to evaluate the polynomials
y = xPowers * coeffs';
end