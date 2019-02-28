function [q, H, Dh, exitflag] = equilibrium_conf_constr(obj, q_0, x, u, w)
%
% q = equilibrium_conf_constr(q_0, x, u, w)
%
% Find the equilibrium configuration given initial guess, tip position of the surface,
% actuation, and external force.
%
% Inputs:
% q_0 is an initial guess of the joint angles.
% x is the desired tip position in R^3.
% u is an actuator input current vector.
% w is an external force vector.
%
% Outputs:
% q is the equilibrium joint angle vector.
% H is the Hessian of the Lagrangian.
% Dh is the Jacobian of the constraint.
% exitflag: 0 means gradient of Lagrangian is smaller than tolerance.
%           1 means step size is smaller than tolerance.
%           -1 means number of iterations exceeds limit.
%           The rest of the flags are from quadprog.
%

    % Stiffness matrix
    K = obj.stiffnessMatrix;
    
    % Initialize variable for optimization algorithm
    q = q_0;
    n = size(q, 1);
    itr = 0;  % Number of iteration
    itr_max = 1000;  % Set number of max iteration
    eps = obj.tolSurfaceConstraint;  % Set tolerance for termination
    B = K;  % Initialize estimated Hessian matrix
    position_constr = @(q)(obj.tip_position(q) - x);  % Tip potition constraints

    while true
        % Store x from previous iteration and update iteration counts
        q_old = q;
        itr = itr + 1;
        
        % Generate the required QP obj.parameterseters for the next iteration
        h = position_constr(q);  % Constraint value at q
        Dh = jacobian_forward_difference(position_constr, q, 1e-8);  % Jacobian of tip constraint
        T = obj.actuator_joint_torque(q, u);  % Actuator torques
        T_d = obj.external_joint_torque(q, w);  % Disturbance torques
        gradf = K*q - T - T_d;  % Gradient of objective function
        
        % Solve the local QP
        options = optimset('Algorithm','interior-point-convex','Display','off');
        [pk, ~, qpflag, ~, lk] = quadprog(B, gradf, [], [], Dh, -h, -pi*ones(n, 1)-q, pi*ones(n,1)-q, zeros(n,1), options);
        
        % Get flag from QP
        if (qpflag == -2 && itr > 1) || qpflag == -3 || qpflag == -6
            exitflag = qpflag;
            break
        end
            
        % Update current points in SQP
        step_size = 0.5;
        q = q + step_size * pk;
        lambda = lk.eqlin;
        
        % Store Gradient of Lagrangian function from previous iteration
        DL = @(q)(K * q + jacobian_forward_difference(position_constr, q, 1e-8)' * lambda - ...
            obj.actuator_joint_torque(q, u) - obj.external_joint_torque(q, w));
        
        % Check for termination
        if norm(DL(q),2) < eps
            exitflag = 1;
            break
        elseif itr == itr_max
            exitflag = 0;
            break
        end
        
        % BFGS update
        B = damped_bfgs(B, q, q_old, DL);
        B = 0.5 * (B + B');  % Force symmetry
    end
    
    % If nothing went wrong, calculate Hessian and constraint gradient at
    % the solution
    if exitflag >= 0
        H = jacobian_forward_difference(DL, q, 1e-4);
        Dh = jacobian_forward_difference(position_constr, q, 1e-8);
    else
        H = [];
        Dh = [];
    end
end
