function [u, jointAngles] = min_contact_(obj, velocity_samples, alpha, q_0, u_0, N_x, x, externalWrenches, frictionCoefficient )

    N_ = size(N_x, 2);

    initialGuess = randn(N_, 1);

    
    if velocity_samples == 0
        lb = -0.1 * ones(size(initialGuess, 1), 1);
        ub = 0.1 * ones(size(initialGuess, 1), 1);
        fun = @(dv)( obj.equilibrium_contact_(q_0, u_0, N_x, dv,externalWrenches, x) );
        constraint = @(dv)point_constraint_(obj,  q_0, u_0, N_x, dv, externalWrenches, x);
        options = optimoptions('fmincon', 'Display','off','TolFun', 1e-06, 'TolCon', 1e-06);
    else
        lb = -0.9 * ones(size(initialGuess, 1), 1);
        ub = 0.9 * ones(size(initialGuess, 1), 1);
        fun =  @(dv)(obj.equilibrium_contact_flow_(velocity_samples, alpha, q_0, u_0, N_x, dv, externalWrenches, x, frictionCoefficient) );
        constraint = @(dv)point_constraint_flow_(obj, velocity_samples, alpha, q_0, u_0, N_x, dv, externalWrenches, x, frictionCoefficient);       
            options = optimoptions('fmincon', 'Display','off','TolFun', 1e-06, 'TolCon', 1e-06);
    end
   
    
    % Use fmincon to minimize potential energy.
    [dv, sigma, exitflag] = fmincon(fun, initialGuess, [], [], [], [], lb, ub, constraint, options);

    u = u_0 + N_x * dv;
    [jointAngles, hessian, lambda, exitflag] = obj.min_potential_energy_conf_const( q_0, u, externalWrenches, [], x, options);
    
end

function [c, ceq] = point_constraint_(obj, q_0, u_0, N_x, dv, externalWrenches, x)
    ceq = [];
%     c = -obj.equilibrium_contact_rh( q_0, u_0, N_x, dv,externalWrenches, x); % sigma has to be positive
    c(1) = -obj.equilibrium_contact_( q_0, u_0, N_x, dv,externalWrenches, x); % for sigma in (0, 0.2)
    c(2) = obj.equilibrium_contact_( q_0, u_0, N_x, dv,externalWrenches, x) - 0.2;
end

function [c, ceq] = point_constraint_flow_(obj, velocity_samples, alpha, q_0, u_0, N_x, dv, externalWrenches, x, frictionCoefficient)
    ceq = [];

    c(1) = -obj.equilibrium_contact_flow_(velocity_samples, alpha, q_0, u_0, N_x, dv, externalWrenches, x, frictionCoefficient); % for P_f, (0, 0.5)
    c(2) = obj.equilibrium_contact_flow_(velocity_samples, alpha, q_0, u_0, N_x, dv, externalWrenches, x, frictionCoefficient) - 0.5;
end

