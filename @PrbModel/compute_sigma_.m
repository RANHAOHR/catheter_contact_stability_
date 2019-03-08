function [sigma_mu, f_c, P_s] = compute_sigma_(obj, velocity_samples, alpha, state, control, disturbances, frictionCoefficient )

    direction_angle = [cos(alpha(1)),cos(alpha(2)),cos(alpha(3))]';

    [F_e] = obj.compute_external_force(velocity_samples, direction_angle, state); %get the external motion caused by the blood flow
    Nsample = size(F_e, 2);
    
    [f_c, sigma_mu, ~] = obj.contact_force_flow_(state, control, disturbances, F_e);
    
    v_safe = sigma_mu(sigma_mu <= frictionCoefficient & sigma_mu >= 0 );
    P_s = size(v_safe,2) / Nsample;

end