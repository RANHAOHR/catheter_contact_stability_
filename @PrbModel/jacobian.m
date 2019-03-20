function [endEffectorJacobian, endEffectorNullspace, quasistaticJacobian, surfaceJacobian] ...
    = jacobian(obj, jointAngles, currents, externalWrenches, hessian)
%
% [endEffectorJacobian, endEffectorNullspace, quasistaticJacobian, surfaceJacobian]
%    = jacobian(obj, jointAngles, currents, externalWrenches, hessian)
%
% Calculate the Jacobian of tip position with respect to input current.
%   J = dp/du = dp/dq * dq/du
% Finite difference is used to calculate the partial derivatives.
% The Jacobian dq/du is calculated according to implicit function theorem.
% 
% Inputs:
% jointAngles is a joint angle vector.
% currents is an actuator current vector.
% externalWrenches is a list of external wrench vector.
% hessian (optional) is the Hessian of the lagrangian of the potential energy
%  minimization problem. It is a part the quasistatic Jacobian, which is based
%  on the implicit function theorem.
%
% Outputs:
% endEffectorJacobian is the Jacobian of the end-effector position with respect
%  to the joint angles.
% endEffectorNullspace is the nullspace of the end-effector Jacobian, where
%  [endEffectorJacobian; endEffectorNullspace] is a full-rank square matrix.
% quasistaticJacobian is the Jacobian of the quasistatic joint angles with
%  respect to the actuation.
% surfaceJacobian is the Jacobian of the spatial mapping (maps from a surface
%  point in R^2 to a point in R^3) of the surface with respect to the 
%  parameterization of the surface.
%
    
    if nargin < 5
        hessian = [];
    end

    % Calculate the end-effector Jacobian analytically.
    endEffectorJacobian = end_effector_jacobian(obj, jointAngles);

    % Calculate the quasi-static Jacobian using the implicit function theorem 
    % and finite difference differentiation.
    quasistaticJacobian = quasistatic_jacobian(obj, jointAngles, currents, ...
        externalWrenches, hessian);
    
    % If the surface is defined and the end-effector is on the surface.
    if (~isempty(obj.surface) && ...
            abs(obj.surface.distance(obj.tip_position(jointAngles))) < ...
            obj.tolSurfaceConstraint)
        endEffectorNullspace = null(endEffectorJacobian)';
        surfaceJacobian = obj.surface.spatial_jacobian(obj.surface.projection(...
            obj.tip_position(jointAngles)));
    else
        endEffectorNullspace = [];
        surfaceJacobian = [];
    end
        

end


function jacobian = end_effector_jacobian(obj, jointAngles)
%
% jacobian = end_effector_jacobian(obj, jointAngles)
%
% This function calculates the Jacobian of the end-effector position with respect
%  to the joint angles.
%
% Input:
% jointAngles is a joint angles vector.
%
% Output:
% jacobian is the Jacobian of the end-effector position with respect to the joint
%  angles.
%

    % Initialize variables.
    jacobian = zeros(3, obj.get_dofs());
    configurations = zeros(4, 4, obj.numJoints + 1);
    configurations(:, :, 1) = eye(4);

    % Calculate the configuration of the joints.
    for i = 1 : 1 : obj.numJoints
        configurations(:, :, i + 1) = configurations(:, :, i) * ...
            se3rot(jointAngles(3*i-2:3*i, 1), [0; 0; obj.jointPositions(1, i)], 1); 
    end

    % Calculate tip position at the current configuration.
    tipPositionHomogeneous = configurations(:, :, end) * [0; 0; obj.length; 1];

    % Calculate the spatial Jacobian.
    for i = 1 : 1 : obj.numJoints
        [rotationAxisX, rotationAxisY, rotationAxisZ] = obj.velocity_axes(jointAngles(3*i-2:3*i));
        adjoint = adj(configurations(:, :, i));
        twistX = adjoint * twistr(rotationAxisX, [0; 0; obj.jointPositions(i)]);
        twistY = adjoint * twistr(rotationAxisY, [0; 0; obj.jointPositions(i)]);
        twistZ = adjoint * twistr(rotationAxisZ, [0; 0; obj.jointPositions(i)]);
        jacobian(:, 3 * i - 2) = [vector2skewsym(twistX(4:6)), twistX(1:3)] * tipPositionHomogeneous;
        jacobian(:, 3 * i - 1) = [vector2skewsym(twistY(4:6)), twistY(1:3)] * tipPositionHomogeneous;
        jacobian(:, 3 * i) = [vector2skewsym(twistZ(4:6)), twistZ(1:3)] * tipPositionHomogeneous;
    end

end

function jacobian = quasistatic_jacobian(obj, jointAngles, currents, externalWrenches, hessian)
%
% jacobian = quasistatic_jacobian(obj, jointAngles)
%
% This function calculates the Jacobian of the quasistatic configuration with
% respect to the currents. The Jacabian is calculated using forward finite 
% different and implicit function theorem.
%
% Inputs:
% jointAngles is a joint angle vector.
% currents is an input current vector.
% externalWrenches is a list of external wrenches.
% epsilon is the surface distance threshold.
%
% Output:
% jacobian is the Jacobian of the quasistatic joint angles with respect to the
%  currents.
%

    % If the Hessian of the Lagrangian of the potential energy minimization 
    % problem is not given, calculate it using the finite difference method.
    if isempty(hessian)
        
        % If the surface is defined and the end-effector is on the surface, use
        % the constrained Lagrangian. Otherwise, use the unconstrained Lagrangian.
        if (~isempty(obj.surface) && abs(obj.surface.distance(obj.tip_position(jointAngles))) < ...
                obj.tolSurfaceConstraint)
            % Assuming that the joint angles is a minimizer of the potential energy,
            % the constraint force is equal to the contact force projected back
            % to the joint space.
            [contactForce, ~] = obj.contact_force(jointAngles, currents, externalWrenches);
            lagrangeMultiplier = norm(contactForce, 2);
            surface_constraint_jacobian = @(q)(obj.surface.distance_jacobian(...
                obj.tip_position(q)) * end_effector_jacobian(obj, q));
            fun = @(q)(obj.stiffnessMatrix * q - ...
                obj.joint_torques(q, currents, externalWrenches) + ...
                lagrangeMultiplier * surface_constraint_jacobian(q)');

        else
            % This is the Lagrangian of the unconstrained catheter.
            fun = @(jointAngles)(obj.stiffnessMatrix * jointAngles - ...
                obj.joint_torques(jointAngles, currents, externalWrenches));
        end
            
        % Calculate the Hessian using the forward difference approach.
        hessian = jacobian_forward_difference(fun, jointAngles, 1e-08);
    end

    % Calculate the quasistatic Jacobian using the implicit function theorem.
    [~, ~, J_u, ~, B, P] = obj.actuation_maps(jointAngles);
    jacobian = -hessian \ (-J_u' * B * P);

end
