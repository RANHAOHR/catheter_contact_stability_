%% Initialize

% Clear workspace.
clear;
close all;
clc;

% Initialize PrbModel
% surface = plane;
numJoints = 4;
parameterFile = './data/parameters_2coilsets.yml';
catheter = PrbModel(parameterFile, numJoints);
% catheter.set_surface(plane); % set the plane with respect to the tip

disturbances = zeros(6, catheter.get_num_joints());

% Calculcate initial configuration
state = [0.5236,0,0,-0.5236,0,0,0.2311,0,0, 0, 0, 0]';
tip_ = catheter.tip_position(state);
control = [0; 0; -1; 0; 0; -1]; % this can be really far away...

theta_vec0 = [[0.5236;0;0], [-0.5236;0;0], [-0.2311;0;0], [0;0;0]];
state = reshape(theta_vec0,12,1);
tip_ = catheter.tip_position(state);
control = [0, 0, -4,0, 0, 3]';

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
% 
% %% Find quasi-static configuration
% 
% % start contact force control
% 
% [J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances);
% 
% J_x = J_e * J_q;
% N_x = null(J_x);
% 
% [control, state] = catheter.min_contact_( state, control, N_x, tip_, disturbances )
% f_c_opt = catheter.contact_force(state, control, disturbances)
% sigma_mu_opt = catheter.compute_contact_ratio(f_c_opt)
% 
% % catheter.tip_position(state)
% % tip_
% catheter.plot_catheter(state, 'red');
% f_c_0
% sigma_mu = catheter.compute_contact_ratio(f_c_0)
% 
% disp('Enter to proceed');
% pause;
%% Add blood flow disturbance
figure(2)
% [velocity_samples] = blood_flow; % get blood flow samples
velocity_samples = [0:0.01:0.9];


beta = pi/2;
direction_angle = [cos(beta),cos(pi/2 - beta),0]';

[F_e] = catheter.compute_external_force(velocity_samples, direction_angle, state);

Nsample = size(F_e, 2);
sigma_mu1 = zeros(1, Nsample);
f_c = zeros(3, Nsample);
for i = 1:Nsample
    [f_c(:,i), ~] = catheter.contact_force_flow_(state, control, disturbances, F_e(:,i));
    sigma_mu1(i) = catheter.compute_contact_ratio(f_c(:,i));
end
figure(2)
plot(velocity_samples, f_c(1,:), 'r-.','LineWidth',2 );
hold on;
plot(velocity_samples, f_c(2,:), 'g-.','LineWidth',2 );
hold on;
plot(velocity_samples, f_c(3,:), 'b-.','LineWidth',2 );
figure(3)
plot(velocity_samples, sigma_mu1, 'y-.','LineWidth',2 );
% 
% v_safe1 = sigma_mu1(sigma_mu1 <= frictionCoefficient & sigma_mu1 >= 0 );
% P_s1 = size(v_safe1,2) / Nsample;
% 
% 
% beta = -pi/2;
% direction_angle = [cos(beta),cos(pi/2 - beta),0]';
% 
% [F_e] = catheter.compute_external_force(velocity_samples, direction_angle, state);
% 
% Nsample = size(F_e, 2);
% sigma_mu2 = zeros(1, Nsample);
% f_c = zeros(3, Nsample);
% for i = 1:Nsample
%     [f_c(:,i), ~] = catheter.contact_force_flow_(state, control, disturbances, F_e(:,i));
%     sigma_mu2(i) = sqrt(f_c(1,i) * f_c(1,i) + f_c(2,i) * f_c(2,i)) / (f_c(3,i));
% end
% 
% v_safe1 = sigma_mu2(sigma_mu2 <= frictionCoefficient & sigma_mu2 >= 0 );
% P_s2 = size(v_safe1,2) / Nsample;