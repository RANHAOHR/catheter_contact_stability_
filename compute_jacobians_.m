function [J_e, J_e_Nullspace, J_q, J_surface]=compute_jacobians_( q_init, l_c, stiffnessMatrix, theta_vec, f_c0_, u, surface_origin, surface_orientation )

    jointAngles = reshape(theta_vec,12,1);
    [J_e, J_e_Nullspace]=compute_endeffector_jacobian_(q_init,jointAngles, l_c);

    [J_q]=compute_quasistatic_jacobian_(l_c, f_c0_, q_init, jointAngles, surface_orientation, stiffnessMatrix, u);
    [J_surface]=surface_jacobian(surface_orientation);
end

function [J_e, J_e_Nullspace]=compute_endeffector_jacobian_(q_init,jointAngles, l_c)

    J_e = zeros(3,12);
    numJoints =4;
    configurations = zeros(4, 4, numJoints + 1);
    configurations(:, :, 1) = eye(4); %the first one does not need to get transformed

    % Calculate the configuration of the joints.
    for i = 1 : 1 : numJoints
        w_vec = jointAngles(3*i-2:3*i, 1); %check if this is a column vector
        g = compute_exponential_map_(w_vec,q_init(:,i)); % x_i
        configurations(:, :, i + 1) = configurations(:, :, i) * g; 
    end

    % Calculate tip position at the current configuration.
    tipPositionHomogeneous = configurations(:, :, end) * [0; 0; l_c; 1];
    
    for i = 1:1:numJoints
        w_vec = jointAngles(3*i-2:3*i, 1);
        [wx, wy, wz] = velocity_axes(w_vec);
        adjoint = adj(configurations(:, :, i));
        
        xi_x = adjoint * [-cross(wx,q_init(:,i));wx];
        xi_y = adjoint * [-cross(wy,q_init(:,i));wy];
        xi_z = adjoint * [-cross(wz,q_init(:,i));wz];
    
        J_e(:, 3*i-2) = [vec_map(xi_x(4:6)), xi_x(1:3)] * tipPositionHomogeneous;
        J_e(:, 3*i-1) = [vec_map(xi_y(4:6)), xi_y(1:3)] * tipPositionHomogeneous;
        J_e(:, 3*i) = [vec_map(xi_z(4:6)), xi_z(1:3)] * tipPositionHomogeneous;

    end
    
    J_e_Nullspace = null(J_e)';
    
end


function [J_q]=compute_quasistatic_jacobian_(l_c, f_c0_, q_init, jointAngles, surface_orientation, stiffnessMatrix, u)
    lagrangeMultiplier = norm(f_c0_, 2);

    surface_constraint_jacobian = @(jointAngles)(distance_jacobian( surface_orientation ) * compute_endeffector_jacobian_(q_init, jointAngles, l_c) ); % get the z axis
    fun = @(q)(stiffnessMatrix * q - ...
                compute_tau(q_init, q, u) + lagrangeMultiplier * surface_constraint_jacobian(q)');

    hessian = jacobian_forward_difference(fun, jointAngles, 1e-08);
    
    [B,P, J_u, ~, ~] = actuation_map(q_init, jointAngles);
    % Calculate the quasistatic Jacobian using the implicit function theorem.
    J_q = -hessian \ (-J_u' * B * P); % equation (10)
end

function [J]=distance_jacobian(surface_orientation)
    J = surface_orientation(:,3)';
end

function [J]=surface_jacobian(surface_orientation)

    J = surface_orientation(:,1:2);
end

