function [jointAngles, hessian, lambda, exitflag] = ...
    min_potential_energy_conf( initialGuess, stiffnessMatrix, gst_0, q_init, l, currents, initialJointAngles, surface_origin, surface_orientation, options)
    % Handle arguments.
    
    % Fill zeros in empty arrays
    if isempty(initialJointAngles) %the configuration at 0
        initialJointAngles = zeros(12, 1);
    end
    
    
    % Set surface constraint.

    constraint = @(jointAngles)surface_constraint(gst_0, q_init, jointAngles, surface_origin, surface_orientation);

    % Define objective function
    fun = @(jointAngles)potential_energy(stiffnessMatrix, q_init, l, jointAngles, currents, initialJointAngles);
    
    % Set lower and upper bounds
    lb = -0.5 * pi * ones(size(initialGuess, 1), 1);
    ub = 0.5 * pi * ones(size(initialGuess, 1), 1);
    
    % Use fmincon to minimize potential energy.
    [jointAngles, ~, exitflag, ~, lambdas, ~, hessian] = fmincon(fun, initialGuess, [], [], [], [], lb, ub, constraint, options);
    
    % Get Lagrangian multiplier for surface constraint
    lambda = lambdas.ineqnonlin;

end

function [c,ceq] = surface_constraint(gst_0, q_init, jointAngles, surface_origin, surface_orientation)
    ceq = [];
    c = (tip_position(gst_0, q_init, jointAngles) - surface_origin)' * surface_orientation(:, 3);  % Tip cannot go through surface
end
