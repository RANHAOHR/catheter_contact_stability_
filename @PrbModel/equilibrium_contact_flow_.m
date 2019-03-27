function [P_f] = equilibrium_contact_flow_(obj,velocity_samples, alpha, q_0, u_0, N_x, dv,externalWrenches, x, frictionCoefficient)
    options = optimoptions('fmincon', 'MaxFunEvals', 1000);

    u = u_0 + N_x * dv;
    [jointAngles, hessian, lambda, exitflag] = obj.min_potential_energy_conf_const(q_0, u, externalWrenches, [], x, options);
    
    % after getting the theta, test contact force
    w_v = [alpha, pi/2 - alpha, pi/2];
    [sigma_mu, f_c, P_s] = obj.compute_sigma_(velocity_samples, w_v, jointAngles, u, externalWrenches, frictionCoefficient );
    P_s
    %If use fmincon here, then use the failure rate instead
    P_f = 1-P_s;
 
end
