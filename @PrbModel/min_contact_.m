function [u, jointAngles] = min_contact_(obj, q_0, u_0, N_x, x, externalWrenches )

    fun =  @(dv)(obj.equilibrium_contact_rh(q_0, u_0, N_x, dv,externalWrenches, x));
    constraint = @(dv)point_constraint(obj,q_0, u_0, N_x, dv,externalWrenches, x);
    N_ = size(N_x, 2);
    options = optimoptions('fmincon', ...
                       'Display', 'off');
    initialGuess = randn(N_, 1);
    lb = -2.5 * ones(size(initialGuess, 1), 1);
    ub = 2.5 * ones(size(initialGuess, 1), 1);
    
    % Use fmincon to minimize potential energy.
    [dv, sigma, exitflag] = fmincon(fun, initialGuess, [], [], [], [], lb, ub, constraint, options);
    exitflag

    u = u_0 + N_x * dv;
    [jointAngles, hessian, lambda, exitflag] = obj.min_potential_energy_conf_const( q_0, u, externalWrenches, [], x, options);
  
%     exitflag
    
end

function [c, ceq] = point_constraint(obj, q_0, u_0, N_x, dv,externalWrenches, x)
    ceq = [];
    c = -obj.equilibrium_contact_rh( q_0, u_0, N_x, dv,externalWrenches, x); % sigma has to be positive
end
