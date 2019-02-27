function [f_c_, N_theta, J_sf, J_C_T_inv] = compute_contact_(gst_0, stiffnessMatrix, q_init, jointAngles, u)

    J_sf =compute_spatial_jacobian_(q_init,jointAngles);

    %compute other forces, gravitational or damping
    N_theta = stiffnessMatrix * jointAngles;
    [tau_u, ~] = compute_tau(q_init, jointAngles, u);
    % In validation, given a reasonable torque, compute f_c and sigma, then add
    % flow, under this torque

    g_tc =  [1,0,0,0;
             0,-1,0,0;
             0,0,-1,0;
             0,0,0,1];

    g = eye(4); %the first one does not need to get transformed

    % Calculate the configuration of the joints.
    for i = 1 : 4
        w_vec = jointAngles(3*i-2:3*i, 1); %check if this is a column vector
        g = g * compute_exponential_map_(w_vec,q_init(:,i));
    end
    g_st = g * gst_0;
    g_sc = g_st * g_tc;

    %compute f_c
    B = [eye(3);zeros(3)];
    J_Cf = B' * adjinv(g_sc) * J_sf;
    J_C_T_inv = pinv(J_Cf');% inv(J_Cf * J_Cf') * J_Cf;

    f_c_ = J_C_T_inv * (tau_u - N_theta);
end