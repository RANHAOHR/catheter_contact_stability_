function [J_CT,force ]= compute_direction(obj, w_v, velocity_samples, state,control, disturbances )
    direction_angle = [cos(w_v(1)),cos(w_v(2)),cos(w_v(3))]';

    [F_e] = obj.compute_external_force(velocity_samples, direction_angle, state); %get the external motion caused by the blood flow
    Nsample = size(F_e, 2);

    [f_c, sigma_mu, jacobian] = obj.contact_force_flow_(state, control, disturbances, F_e);
    J_CT = pinv(jacobian');
    
    K_theta = -obj.stiffnessMatrix * state;
    K_theta = repmat(K_theta, [1,Nsample]);
    tau = obj.joint_torques(state, control, disturbances);
    tau = repmat(tau, [1,Nsample]);

    force = K_theta + tau + F_e; % should be a 12*n matrix with forces in columns
end