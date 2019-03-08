%% Initialize

% Clear workspace.
clear;
close all;
clc;
set(0,'defaultfigurecolor','w');
% Initialize PrbModel
% surface = plane;
numJoints = 4;

parameterFile = './data/parameters_2coilsets.yml';
catheter = PrbModel(parameterFile, numJoints);
% catheter.set_surface(plane); % set the plane with respect to the tip

disturbances = zeros(6, catheter.get_num_joints());

% Calculcate initial configurations
% state = [0.5236,0,0,-0.5236,0,0,0.2311,0,0, 0, 0, 0]';
% tip_ = catheter.tip_position(state);
% control = [0; 0; 3; 0; 0; -3]; % this can be really far away...

state = [0.0725, -0.0514, 1.1160, 0.1259, -0.2842, 1.1149, -0.2702, -0.1743, -1.5707, -0.2305, 0.2968, 0.0]';
tip_ = catheter.tip_position(state);
control = [0.6544, -2.0530, -1.5000, -0.08723, 1.5228, -1.9458]'; % this can be really far away...

% state = [-0.5236,0,0,0.5236,0,0,-0.2311,0,0, 0, 0, 0]';
% tip_ = catheter.tip_position(state);
% control = [0.0; 0; 1; 0.0; -0.1; 1.5]; % this can be really far away...

% theta_vec0 = [[0.5236;0;0], [-0.5236;0;0], [-0.2311;0;0], [0;0;0]];
% state = reshape(theta_vec0,12,1);
% tip_ = catheter.tip_position(state);
% control = [0, 0, -4,0, 0, 3]';

% Create a plane
origin = [0; 0; tip_(3)]; %% just to gurantee on the surface plane now...
orientation = so3rot([1; 0; 0], pi);
frictionCoefficient = 0.2;

plane = Plane(origin, orientation, frictionCoefficient);
catheter.set_surface(plane);

% catheter.plot_catheter(state, 'blue');
% options = optimoptions('fmincon', 'MaxFunEvals', 1000);
% [state, hessian, lambda, exitflag] = catheter.min_potential_energy_conf_const(...
%             state, control, disturbances, [], tip_, options)
% tip_ = catheter.tip_position(state)
% catheter.plot_catheter(state, 'red');

% Print results
disp('contact_force = ');
f_c_0 = catheter.contact_force(state, control, disturbances)
catheter.tip_position(state)
catheter.plot_catheter(state, 'blue');
pause;


%% Find quasi-static configuration

% start contact force control

[J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances); %endEffectorNullspace here is the transpose fo the null space...

J_x = J_e * J_q;
N_x = null(J_x);

[velocity_samples] = blood_flow;
% velocity_samples = 0;

alpha_ = 4.01;
alpha = [alpha_, pi/2 - alpha_, pi/2];

[sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );
disp('initial:');
P_s

% [control, state] = catheter.min_contact_( state, control, N_x, tip_, disturbances )
% f_c_opt = catheter.contact_force(state, control, disturbances)
% sigma_mu_opt = catheter.compute_contact_ratio(f_c_opt)

[control, state] = catheter.min_contact_(velocity_samples, alpha, state, control, N_x, tip_, disturbances, frictionCoefficient );
[sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );

% catheter.tip_position(state)
% tip_
catheter.plot_catheter(state, 'red');
f_c_0
sigma_mu = catheter.compute_contact_ratio(f_c_0)

disp('Enter to proceed');
pause;

%% evaluate blood flow disturbance
% % % 
% alpha_range = [0, 2*pi]';
% w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
% catheter.velocity_angle_analysis(alpha_range, w_v,  state, control, disturbances, frictionCoefficient);

% velocity_samples = 0.1;
% alpha = 0;
% w_v_alpha = [alpha, pi/2 - alpha, pi/2]';
% externalWrenches = zeros(6,4);
% [F_e] = catheter.compute_external_force(velocity_samples, w_v_alpha, state);
% [f_c, sigma_mu, jacobian] = catheter.contact_force_flow_(state, control, externalWrenches, F_e)