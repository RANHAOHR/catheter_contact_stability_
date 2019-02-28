function [jointAngles, hessian, lambda, exitflag] = ...
    min_potential_energy_conf_const(obj, initialGuess, currents, externalWrenches, initialJointAngles, x, options)
%
% [jointAngles, hessian, lambda, exitflag] =
%   min_potential_energy(initialGuess, currents, externalWrenches, initialJointAngles, options)
%
% Calculate equilibrium configuration.
%
% Inputs:
% initialGuess is an initial guess of the joint angle vector.
% currents is an actuator input current vector.
% externalWrenches is an array of external wrenches.
% initialJointAngles is a joint angle vector at the initial configuration.
% options is Matlab optimization option struct.
%
% Outputs:
% jointAngles is the equilibrium joint angle vector.
% hessian is the Hessian of the Lagrangian.
% lambda is the Lagrange multipliers.
% exitflag: 0 means gradient of Lagrangian is smaller than tolerance
%           1 means step size is smaller than tolerance
%           -1 means number of iterations exceeds limit
%           The rest of the flags are from quadprog
%

    % Handle arguments.
    if nargin < 3
        error('Not enough input arguments.');
    end
    
    if nargin < 6
        options = optimoptions('fmincon', ...
                               'Display', 'off', ...
                               'TolFun', obj.tolPotentialEnergy, ...
                               'TolCon', obj.tolSurfaceConstraint);
    end
    
    if nargin < 5
        initialJointAngles = [];
    end
    
    if nargin < 4
        externalWrenches = [];
    end
    
    % Fill zeros in empty arrays
    if isempty(initialJointAngles)
        initialJointAngles = zeros(obj.get_dofs(), 1);
    end
    
    if isempty(externalWrenches)
        externalWrenches = zeros(6, obj.numJoints, 1);  % Potential energy minimization currently does not consider external wrenches.
    end
    
    % Set surface constraint.
    if isempty(obj.surface)
        constraint = [];
    else
        constraint = @(jointAngles)surface_constraint(obj, jointAngles,x);
    end

    % Define objective function
    fun = @(jointAngles)obj.potential_energy(jointAngles, currents, externalWrenches, initialJointAngles);
    
    % Set lower and upper bounds
    lb = -0.5 * pi * ones(size(initialGuess, 1), 1);
    ub = 0.5 * pi * ones(size(initialGuess, 1), 1);
    
    % Use fmincon to minimize potential energy.
    [jointAngles, ~, exitflag, ~, lambdas, ~, hessian] = fmincon(fun, initialGuess, [], [], [], [], lb, ub, constraint, options);
    
    % Get Lagrangian multiplier for surface constraint
    lambda = lambdas.ineqnonlin;

end

function [c, ceq] = surface_constraint(obj, jointAngles, x)
    ceq = obj.tip_position(jointAngles) - x;
    c = [];  % Tip cannot go through surface
end
