function [sigma_mu, f_c, P_s] = compute_sigma_(obj, velocity_samples, beta, state, control, disturbances, frictionCoefficient )

    direction_angle = [cos(beta(1)),cos(beta(2)),cos(beta(3))]';

    [F_e] = obj.compute_external_force(velocity_samples, direction_angle, state);

    Nsample = size(F_e, 2);
    sigma_mu = zeros(1, Nsample);
    f_c = zeros(3, Nsample);
    for i = 1:Nsample
        [f_c(:,i), ~] = obj.contact_force_flow_(state, control, disturbances, F_e(:,i));
        sigma_mu(i) = obj.compute_contact_ratio(f_c(:,i));
    end

    v_safe = sigma_mu(sigma_mu <= frictionCoefficient & sigma_mu >= 0 );
    P_s = size(v_safe,2) / Nsample;

end