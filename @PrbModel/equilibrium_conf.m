function [q, H, exitflag] = equilibrium_conf(obj, initialGuess, currents, externalWrenches, options)
%
% [q, H, exitflag] = equilibrium_conf(q0, u, w)
%
% Find the equilibrium configuration given initial guess, actuation,
% and external force.
%
% Inputs:
% initialGuess is an initial guess of the joint angle vector.
% currents is an actuator input current vector.
% externalWrenches is an array of external wrenches.
% options is Matlab optimization option struct.
%
% Outputs:
% q is the joint angles.
% H is the Hessian of the Lagrangian.
% exitflag: 0 means gradient of Lagrangian is smaller than tolerance.
%           1 means step size is smaller than tolerance.
%           -1 means number of iterations exceeds limit.
%           The rest of the flags are from quadprog.
%    

    % Optimization options is not given
    if nargin == 4
        options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off');
    end

    % obj.parameterseters
    K = obj.stiffnessMatrix;
    % Set the initial guess
    q = initialGuess;
    % DOF
    n = size(q, 1);
    % Number of iteration
    itr = 0;
    % Set number of max iteration
    itr_max = 1000;
    % Set tolerance for the gradient of the Lagrangian
    eps = obj.tolEquilibrium;
    % Generate the required QP obj.parameterseters
    % Initialize BFGS matrix
    B = K;
    % Step size
    step_size = 1e0;

    % Iterate until one of the termination criteria is met
    while true
        % Store previous solution and update iteration counts
        q_old = q;
        itr = itr + 1;
        
        % Calculate gradient for QP subproblem
        T = obj.actuator_joint_torque(q, currents);  % Actuator torques
        T_d = obj.external_joint_torque(q, externalWrenches);  % Disturbance torques
        gradf = K * q - T - T_d;
        
        % Solve the local QP
        [pk, ~, qpflag] = quadprog(B, gradf, [], [], [], [], -0.5 * pi * ones(n, 1) - q, 0.5 *pi * ones(n, 1) - q, zeros(n, 1), options);
        
        % Get flag from QP
        if (qpflag == -2 && itr > 1) || qpflag == -3 || qpflag == -6
            exitflag = qpflag;
            break
        end
        
        % Update current point in SQP
        q = q + step_size * pk;  % Line search would be nice here
        
        % Store Gradient of Lagrangian function from previous iteration 
        DL = @(q)(K * q - obj.actuator_joint_torque(q, currents) - obj.external_joint_torque(q, externalWrenches));
        
        % Check for termination
        if norm(DL(q),2) < eps
            exitflag = 1;
            break
        elseif itr == itr_max
            exitflag = 0;
            break
        end
        
        % Perform damped BFGS update
        B = damped_bfgs(B, q, q_old, DL);
        B = 0.5 * (B + B');  % Force symmetry
        
    end
    
    % 
    if exitflag >= 0
        H = jacobian_forward_difference(DL, q, 1e-4);
    else
        H = [];
    end
    
end
