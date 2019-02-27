function [sigma_mu, ps, P_s1, P_s2, P_s3, P_s4]=catheter_kinematics(direction_angle, velocity_samples, theta_vec, f_c0)
    l_c = 0.104; %m
    l = [0.01655,0.2181,0.2893,0.2455,0.2306] * l_c;

    % start Jacobian computation 
    q4 = [0, 0, sum(l(1:4))]';
    q3 = [0, 0, sum(l(1:3))]';
    q2 = [0, 0, sum(l(1:2))]';
    q1 = [0, 0, l(1)]';

    theta_vec0 = [[0.5236;0;0], [-0.5236;0;0], [-0.2311;0;0], [0;0;0]];
    jointAngles0 = reshape(theta_vec0,12,1);
    u_0 = [0, 0, -9 * 10^(-3),0, 0, 8 * 10^(-3)]';
    numJoints = 4;
    outerRadius = 1.5875 * 10^(-3);
    innerRadius = 0.9906 * 10^(-3);
    youngsModulus= 5.3948 * 10^(6);
    shearModulus= 2.3881 * 10^(6);
    moment_inertia = 0.25 * pi * (outerRadius^4 - innerRadius^4);
    % Polar moment of inertia of area
    polar_moment_inertia = 0.5 * pi * (outerRadius^4 - innerRadius^4);
    % Stiffness vector
    stiffnesses = zeros(1, 2 * numJoints);
    
    for i = 1 : numJoints
        stiffnesses(1, 2 * i - 1) = youngsModulus * moment_inertia / l(i + 1);
        stiffnesses(1, 2 * i) = shearModulus * polar_moment_inertia / l(i + 1);
    end

    stiffnessMatrix = diag(kron(eye(4), [1, 0; 1, 0; 0, 1]) * stiffnesses');

    gst_0 = [1,0,0,0;
             0,1,0,0;
             0,0,1,l_c;
             0,0,0,1];

    q_init = [q1, q2, q3, q4 ];
 
    exp_cell_ = cell(1, 4);
    R_cell_ = cell(1,4);
    for i = 1:4
        w_vec = theta_vec(:,i); %check if this is a column vector
        exp_cell_(i) = {compute_exponential_map_(w_vec,q_init(:,i))}; % x_i
        g = cell2mat(exp_cell_(i));
        R_cell_(i) = {g(1:3,1:3)};
    end

    g_st = eye(4);
    for i =1:4
        g_st = g_st * cell2mat(exp_cell_(i));
    end

    g_st = g_st * gst_0;
%     J_b_st=compute_body_jacobian_(gst_0, q_init, exp_cell_, theta_vec, 4);
    
    J_sf =compute_spatial_jacobian_(q_init,jointAngles0);

    %compute other forces, gravitational or damping
    N_theta = stiffnessMatrix * jointAngles0;
    [tau_u, ~] = compute_tau(q_init, jointAngles0, u_0)
    % In validation, given a reasonable torque, compute f_c and sigma, then add
    % flow, under this torque

    g_tc =  [1,0,0,0;
             0,-1,0,0;
             0,0,-1,0;
             0,0,0,1];

    g_sc = g_st * g_tc;

    %compute f_c
    B = [eye(3);zeros(3)];
    J_Cf = B' * adjinv(g_sc) * J_sf;
    J_C_T_inv = pinv(J_Cf');% inv(J_Cf * J_Cf') * J_Cf;

    f_c0_ = J_C_T_inv * (tau_u - N_theta);
    
      
    surface_origin = tip_position(gst_0, q_init, jointAngles0);  
    surface_orientation = [1,0,0;
             0,-1,0;
             0,0,-1];
    %%%%%% start compute the current...
%     v_N = J_e_Nullspace * delta_theta;
%     step = 1000;
%     dv_N = v_N / step;

    options = optimoptions('fmincon', 'MaxFunEvals', 1000);
    
    dx = [0.1; 0.0];
    numSteps = 3;
    states = zeros(12, numSteps+1);
    states(:,1) = jointAngles0;
    u_ = u_0;
    u_vec = u_0;
    initialJointAngles = zeros(12, 1);
    figure(7);
    for j = 1:numSteps
        [J_e, J_e_Nullspace, J_q, J_surface]=compute_jacobians_( q_init, l_c, stiffnessMatrix, states(:,j) , f_c0_, u_, surface_origin, surface_orientation );

        dz = pinv([J_q, -J_e_Nullspace']) * pinv(J_e) * J_surface * dx;
        du = dz(1:6);       
        u_ = u_ + du;
        
        [states(:,j+1), hessian, lambda, exitflag] = min_potential_energy_conf( states(:,j), stiffnessMatrix,gst_0, q_init, l, u_, initialJointAngles, surface_origin, surface_orientation, options)
        fw_kinematics(l_c, l, states(:,j+1), 'k');
        pause;
        u_vec = [u_vec,u_];
    end
%     N_theta = zeros(12,1);
%     tau_u = N_theta + J_Cf' * f_c0


    %compute blood flow drag
    [F_e]=compute_external_force(velocity_samples, direction_angle, exp_cell_, R_cell_, l, g_st);

    Nsample = size(F_e,2);
    sigma_mu = zeros(1, Nsample);
    f_c = zeros(3, Nsample);
    for i = 1:Nsample
        f_c(:,i) = J_C_T_inv * (J_sf' * F_e(:,i) + tau_u - N_theta);
        sigma_mu(i) = sqrt(f_c(1,i) * f_c(1,i) + f_c(2,i) * f_c(2,i)) / (f_c(3,i));
    end
sigma_mu_ = sqrt(f_c0_(1) * f_c0_(1) + f_c0_(2) * f_c0_(2)) / f_c0_(3)
%     subplot(2,1,2);
%     for i = 1: Nsample
%         line([sigma_mu(i),sigma_mu(i)],[0,0.5]);
%         hold on;
%     end
%     set(gca,'ytick',[]) %need not show the y axis.
%     axis([-1 1 0 1]);
%     xlabel('Distribution of $\sigma_\mu$ samples m/s','Interpreter','latex');
% Setfontsize = 15;
% set(get(gca,'Xlabel'),'FontSize',Setfontsize);
% set(get(gca,'Ylabel'),'FontSize',Setfontsize);    
    mu_s = 0.2;

    v_safe1 = sigma_mu(sigma_mu <= mu_s & sigma_mu >= 0 );
    P_s1 = size(v_safe1,2);
    
    v_safe3 = sigma_mu(sigma_mu > mu_s & sigma_mu <= 0.2 );
    P_s2 = size(v_safe3,2);

    v_safe4 = sigma_mu(sigma_mu > 0.2 );
    P_s3 = size(v_safe4,2);
    
    [~,ind] = find(sigma_mu < 0);
    P_s4 = velocity_samples(ind);
    
    ps = P_s1/ Nsample;


end


