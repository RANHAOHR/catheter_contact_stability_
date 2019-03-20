function [J_cu, J_ctheta, J_cq, J_tz] = compute_contact_jacbobian(obj, jointAngles, currents, Fe, J_e, disturbances)

    %calculate it using the finite difference method.       
    % If the surface is defined and the end-effector is on the surface, use
    % the constrained Lagrangian. Otherwise, use the unconstrained Lagrangian.
    tip_ = obj.tip_position(jointAngles)
    if (~isempty(obj.surface) && abs(obj.surface.distance(obj.tip_position(jointAngles))) < 100)
        % Assuming that the joint angles is a minimizer of the potential energy,
        % the constraint force is equal to the contact force projected back
        % to the joint space.
        N_k = null(J_e);
        dv = zeros(1, size(N_k, 2))';
        
        fun_u = @(u)( pinv(contact_jacobian(obj, jointAngles)') * (obj.joint_torques(jointAngles, u, disturbances) + Fe - obj.stiffnessMatrix * jointAngles) );
        fun_theta = @(q)( pinv(contact_jacobian(obj, q)') * (obj.joint_torques(q, currents, disturbances) + Fe - obj.stiffnessMatrix * q) );
        fun_z = @(z)( pinv(contact_jacobian(obj, N_k * z + jointAngles)') * (obj.joint_torques(N_k * z + jointAngles, currents, disturbances) + Fe - obj.stiffnessMatrix *( N_k * z + jointAngles) ) );
    
        fun_q = @(q)( obj.stiffnessMatrix * q - ...
                obj.joint_torques(q, currents, disturbances) ); 
        

%         fun_kz = @(z)( obj.stiffnessMatrix * (N_k * z + jointAngles) - ...
%                 obj.joint_torques(N_k * z + jointAngles, currents, disturbances) ); 
%         fun_ku = @(u)( obj.stiffnessMatrix * (N_k * dv + jointAngles) - ...
%                 obj.joint_torques(N_k * dv + jointAngles, u, disturbances) ); 
    else
        error('No contact detetcted.');
        return
    end
%     J_ku = jacobian_forward_difference(fun_ku, currents, 1e-08);
%     J_kz = jacobian_forward_difference(fun_kz, dv, 1e-08);
%     
%     J_cq = -J_kz \ J_ku;
% %     fun_z = @(z)( N_k * z + jointAngles); 
%     J_tz = jacobian_forward_difference(fun_z, dv, 1e-08);
J_tz = 0;
    %using the forward difference approach.
    J_cu = jacobian_forward_difference(fun_u, currents, 1e-08);
    J_ctheta = jacobian_forward_difference(fun_theta, jointAngles, 1e-08);
    hessian_q = jacobian_forward_difference(fun_q, jointAngles, 1e-08);
   
    lambda = ( hessian_q ) / J_e;
    fun_ = @(q)(obj.stiffnessMatrix * q - ...
                obj.joint_torques(q, currents, disturbances) + lambda * (tip_ - obj.tip_position(q)) ); 
            
    hessian = jacobian_forward_difference(fun_, jointAngles, 1e-08);
    % Calculate the quasistatic Jacobian using the implicit function theorem.
    [~, ~, J_u, ~, B, P] = obj.actuation_maps(jointAngles);
    J_cq = -hessian_q \ (-J_u' * B * P);
    
   
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
    
    jacobian = null(jacobian);
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

function [f_c] = contact_force_(obj, jointAngles, currents, externalWrenches, F_e)
    Nsample = size(F_e, 2);
    [jacobian] = contact_jacobian(obj, jointAngles);
    K_theta = -obj.stiffnessMatrix * jointAngles;
    K_theta = repmat(K_theta, [1,Nsample]);
    tau = obj.joint_torques(jointAngles, currents, externalWrenches);
    tau = repmat(tau, [1,Nsample]);

    force_ = K_theta + tau + F_e; % should be a 12*n matrix with forces in columns

    f_c = pinv(jacobian') * force_; %the contact force'
end