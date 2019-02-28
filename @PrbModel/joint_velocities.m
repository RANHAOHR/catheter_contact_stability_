function velocities = joint_velocities(obj, angles, currents, externalWrenches)
%
% joint_velocity = joint_velocities(angles, currents, externalWrenches)
%
% This function calculates joint angle velocities given joint angles,
% currents, and external forces. The calculation of Lagrange multiplier (lambda)
% is based on Equation 6.6 in [1].
%
% Inputs:
% angles is a joint angle vector.
% currents is a current vector.
% externalWrenches is an external wrench vector.
%
% Output:
% velocities is a joint angles velocity vector.
%

    jointTorques = obj.joint_torques(angles, currents, externalWrenches);
    dampingInv = inv(obj.damping_matrix(angles));

    % If the surface is defined and the end-effector is on the surface.
    if (~isempty(obj.surface) && abs(obj.surface.distance(obj.tip_position(angles))) < obj.tolSurfaceConstraint)
        
        % Calculate contact force.
        [contactForce, contactJacobian] = obj.contact_force(angles, currents, externalWrenches);
        
        % Check if contact force is in friction cone or not
        if sqrt(contactForce(1)^2 + contactForce(2)^2) <= ...
                obj.surface.get_friction_coefficient() * contactForce(3)
            frictionForce = [-contactForce(1); -contactForce(2); 0];
        else
            % Calculate tip velocity.
            tipVelocity = contactJacobian * dampingInv * (-obj.stiffnessMatrix * ...
                angles + jointTorques);
            tipVelocity(3) = 0;
            
            % Friction force
            if norm(tipVelocity, 2) < 1e-16
                frictionForce = zeros(3, 1);
            else
                frictionForce = -obj.surface.get_friction_coefficient() * ...
                    contactForce(3) * tipVelocity / norm(tipVelocity, 2);
            end
        end
        
        % Calculate joint torques from friction.
        frictionJointTorques = contactJacobian' * frictionForce;
        
        % Calculate the Lagrange multiplier.
        constraint = @(q)obj.surface.distance(obj.tip_position(q));
        constraintJacobian = jacobian_forward_difference(constraint, angles, 1e-8);
        lambda = (constraintJacobian * dampingInv * constraintJacobian') \ ...
            (constraintJacobian * dampingInv * (-obj.stiffnessMatrix * angles + ...
            jointTorques + frictionJointTorques));
        
        % Calculate the difference between internal and external joint torques.
        jointTorquesDiff = -obj.stiffnessMatrix * angles - lambda * constraintJacobian' ...
            + jointTorques + frictionJointTorques;
    else
        % Calculate the difference between internal and external joint torques.
        jointTorquesDiff = -obj.stiffnessMatrix * angles + jointTorques;
    end
    
    velocities = dampingInv * jointTorquesDiff;

end
