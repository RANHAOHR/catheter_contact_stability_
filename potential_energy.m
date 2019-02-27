function f = potential_energy(stiffnessMatrix, q_init, l, jointAngles, u, initialJointAngles)
    
    gravity = [0,0,9.8]';
    magneticField = [0,3,0]';
    coilMasses = 10^(-2);
    % Initialize objective function value with internal stiffness potential
    f = 0.5 * (jointAngles - initialJointAngles)' * stiffnessMatrix * (jointAngles - initialJointAngles);
    
    % Add potential energy from coils
    actuatorConfigurations = zeros(4,4,2);
   
    w_vec = jointAngles(1:3, 1);
    g = compute_exponential_map_(w_vec, q_init(:,1));   
    l_a_z = [0.044,0.0752]';
    
    configuration = eye(4);
    for i = 1:2
        configuration(3,4) = l_a_z(i);
        w_vec = jointAngles(3*(i+1)-2:3*(i+1), 1);
        g = g * compute_exponential_map_(w_vec, q_init(:,(i+1))) ;   
        actuatorConfigurations(:,:,i) = g * configuration; %actuator1 = e1 e2 laz_1, actuator2 = e1 e2 e3 laz_2
        
    end
    turnAreaMatrices = 1e-3 * [30 * 46.6650, 30 * 46.6650, 100 * 11.251791] .* eye(3);
    for i = 1 : 1 : 2
        % Magnetic moment aligned with spatial frame
        magneticMoment = actuatorConfigurations(1:3, 1:3, i) * turnAreaMatrices * ...
            u(3 * (i - 1) + 1 : 3 * (i - 1) + 3);
        % Update potential enerfy
        f = f - magneticMoment' * magneticField - ...
            coilMasses * gravity' * actuatorConfigurations(1:3, 4, i);
    end
    
    % Add potential energy from link weights
    link_origin_0 = zeros(3,5);
    link_origin_0(3,1) = 0.5 * l(1);
    for i = 2:5
        link_origin_0(3,i) = 0.5 * l(i) + sum(l(1:i-1));
    end
    link_origin = zeros(4,5);
    g = eye(4);
    for i =1:5      
        point_vec = [squeeze(link_origin_0(:,i));1];
        link_origin(:,i) = g * point_vec; % the link origins in spatial frame
        if i < 5
            w_vec = jointAngles(3*(i)-2:3*(i), 1);
            g = g * compute_exponential_map_(w_vec, q_init(:,i)) ;    
        end
    end
    mass_vec = [0.01655,0.2181,0.2893,0.2455,0.2306] * 1.2 * 10^(-3); %kg
    for i = 1 : 4
        % Update potential energy
        f = f - mass_vec(1, i) * gravity' * link_origin(1:3,i+1);
    end
end