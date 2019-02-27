function [J_b, J_u]=compute_actuator_jacobian_(g_0,q_init, jointAngles, link_ind)

    J_b = zeros(6,12);
    
    g = g_0;
    for i = link_ind:-1:1
        w_vec = jointAngles(3*i-2:3*i, 1);
        g = compute_exponential_map_(w_vec, q_init(:,i)) * g;
        [wx, wy, wz] = velocity_axes(w_vec);

        xi_x = [-cross(wx,q_init(:,i));wx];
        xi_y = [-cross(wy,q_init(:,i));wy];
        xi_z = [-cross(wz,q_init(:,i));wz];
    
        J_b(:, 3*i-2) = adjinv(g) * xi_x;
        J_b(:, 3*i-1) = adjinv(g) * xi_y;
        J_b(:, 3*i) = adjinv(g) * xi_z;

    end
    
    J_u = J_b(4:6,:); % for later computing torque, with corresponding B
end


