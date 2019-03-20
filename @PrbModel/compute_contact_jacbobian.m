function [J_cu, J_ctheta, J_cq] = compute_contact_jacbobian(obj, jointAngles, currents, Fe, J_e, f_c, disturbances)

    %calculate it using the finite difference method.       
    % If the surface is defined and the end-effector is on the surface, use
    % the constrained Lagrangian. Otherwise, use the unconstrained Lagrangian.
    tip_ = obj.tip_position(jointAngles);
    if (~isempty(obj.surface) && abs(obj.surface.distance(obj.tip_position(jointAngles))) < 100)
        % Assuming that the joint angles is a minimizer of the potential energy,
        % the constraint force is equal to the contact force projected back
        % to the joint space.
        
        fun_u = @(u)( pinv(contact_jacobian(obj, jointAngles)') * (obj.joint_torques(jointAngles, u, disturbances) + Fe - obj.stiffnessMatrix * jointAngles) );
        fun_theta = @(q)( pinv(contact_jacobian(obj, q)') * (obj.joint_torques(q, currents, disturbances) + Fe - obj.stiffnessMatrix * q) );

%         obj.surface.distance_jacobian( obj.tip_position(jointAngles))
%         surface_constraint_jacobian = (obj.surface.distance_jacobian(...
%                 obj.tip_position(jointAngles)) * end_effector_jacobian(obj, jointAngles))
%             end_effector_jacobian(obj, jointAngles)
%             1/ norm(f_c) *
        fun_q = @(q)( obj.stiffnessMatrix * q - ...
                obj.joint_torques(q, currents, disturbances) + 2 / norm(f_c) * (obj.surface.distance_jacobian(...
                obj.tip_position(q)) * end_effector_jacobian(obj, q))' );
            
  
    else
        error('No contact detetcted.');
        return
    end

    %using the forward difference approach.
    J_cu = jacobian_forward_difference(fun_u, currents, 1e-08);
    J_ctheta = jacobian_forward_difference(fun_theta, jointAngles, 1e-08);
    hessian_q = jacobian_forward_difference(fun_q, jointAngles, 1e-08);
    

    
    
   
    % Calculate the quasistatic Jacobian using the implicit function theorem.
    [~, ~, J_u, ~, B, P] = obj.actuation_maps(jointAngles);
    J_cq = -hessian_q \ (-J_u' * B * P);
    
   
end

function [J] = contact_jacobian(obj, q)
%
% J = contact_jacobian(q, param)
%
% Calculate the contact Jacobian (generally know as hand Jacobian).
% See Euation 5.14 in [1] for more details.
%
% Input:
% q is the joint angle vector.
%
% Output:
% J is the Hand Jacobian ().
%
    % Get parameters.
    g_sc = obj.surface.tangent_frame(obj.surface.projection(obj.tip_position(q))); % Assume tip is on the surface
    
    g_sf = eye(4);
    p0 = obj.jointPositions;

    % Define wrench basis.
    B_c = [1 0 0; 0 1 0; 0 0 1; zeros(3)];

    % Initialize Jacobian from base to tip frame.
    J_sf = zeros(6, 3 * obj.numJoints);

    % Calculate the finger Jacobian.
    for n = 1:1:obj.numJoints
        w = q(3 * n - 2 : 3 * n);
        p = [0; 0; p0(1, n)];
        g = se3rot(w, p, 1);

        % Calculate the twists analytically.
        [w1, w2, w3] = obj.velocity_axes(w);
        xi1 = twistr(w1, p);
        xi2 = twistr(w2, p);
        xi3 = twistr(w3, p);

        % Calculate the finer Jacobian.
        J_sf(:, 3 * n - 2) = adj(g_sf) * xi1;
        J_sf(:, 3 * n - 1) = adj(g_sf) * xi2;
        J_sf(:, 3 * n) = adj(g_sf)* xi3;
        g_sf = g_sf * g;
    end

    % Assemble the contact (or hand) Jacobian eq5.14 MLS
    J = B_c' * adjinv(g_sc) * J_sf;

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