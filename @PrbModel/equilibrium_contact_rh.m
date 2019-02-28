function [sigma_mu, jointAngles] = equilibrium_contact_rh(obj, q_0, u_0, N_x, dv,externalWrenches, x)
    options = optimoptions('fmincon', 'MaxFunEvals', 1000);

    u = u_0 + N_x * dv;
    [jointAngles, hessian, lambda, exitflag] = obj.min_potential_energy_conf_const(q_0, u, externalWrenches, [], x, options);
    % after getting the theta, test contact force
    [f_c_, ~] = obj.contact_force(jointAngles, u, externalWrenches);
    f_c_
    sigma_mu = sqrt(f_c_(1) * f_c_(1) + f_c_(2) * f_c_(2)) / f_c_(3)
    
end
