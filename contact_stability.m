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

% state = [-0.0374;0.04693;0.9673;-0.05602;0.3445;0.96563;0.02618;0.08565;0.6394;-0.3222;-0.2276;-0.00058];
% tip_ = catheter.tip_position(state);
% control = [0.1427;-0.1476;1.0;0.0292;-0.3845;1.4167]; % this can be really far away...

% state = [-0.04768;0.03084;0.9315;-0.08019;0.37607;0.93007;-0.021768;0.028814;0.68311;-0.29687;-0.22187;0.0003];
% tip_ = catheter.tip_position(state);
% control = [0.1710;-0.1237;1.2762;0.03464;-0.3939;1.4096]; % this can be really far away...

% state = [-0.0737020339466517;-0.0218040664877389;0.557950048047828;-0.255235186705955;0.337803859876394;0.557141793151812;-0.179753774914578;-0.0895241380738602;1.55912523055963;-0.240881919228409;-0.175779886865597;1.43491696547212e-05];
% tip_ = catheter.tip_position(state);
% control = [-0.318923394305813;0.247555477222941;1.35576833186285;0.0412604294881832;-1.31774443872355;1.48556110031736]; % this can be really far away...

% Create a plane
origin = [0; 0; tip_(3)]; %% just to gurantee on the surface plane now...
orientation = so3rot([1; 0; 0], pi);
frictionCoefficient = 0.2;

plane = Plane(origin, orientation, frictionCoefficient);
catheter.set_surface(plane);

% Print results
figure(1);
disp('contact_force = ');
state_0 = state;
f_c_0 = catheter.contact_force(state, control, disturbances)
sigma_mu_0 = catheter.compute_contact_ratio(f_c_0)
catheter.tip_position(state_0)
catheter.plot_catheter(state_0, 'blue');
pause;

%% start contact force searching

% [J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances); %endEffectorNullspace here is the transpose fo the null space...
% J_x = J_e * J_q;
% N_x = null(J_x);
% 
% velocity_samples = 0;
% alpha = 0;
% 
% [control, state] = catheter.min_contact_(velocity_samples, alpha, state, control, N_x, tip_, disturbances, frictionCoefficient );
% [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );

%% start contact force control 
% figure(2);
% [velocity_samples] = blood_flow;
% 
% alpha_ = pi/2;
% alpha = [alpha_, pi/2 - alpha_, pi/2];
% 
% [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );
% pause;
% 
% [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );

% alpha = [3.6652:0.1:4.3633]';
% disp('Initial P_s:');
% P_s = 0; %consider all P_s of all alphas
% for i = 1: size(alpha,1)
%     w_v = [alpha(i), pi/2 - alpha(i), pi/2];
%     [sigma_mu, f_c, P_s_i] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
%     P_s = P_s + P_s_i;
% end
% P_s = P_s / size(alpha,1)
% 
% alpha_range = [0, 2*pi]';
% w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
% catheter.plot_P_s(alpha_range, w_v, state, control, disturbances, frictionCoefficient)
% pause;

% figure(3);
% Ntrial = 5;
% for k = 1:Ntrial
%     tip_ = catheter.tip_position(state);
%     [J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances); %endEffectorNullspace here is the transpose fo the null space...
% 
%     J_x = J_e * J_q;
%     N_x = null(J_x);
% 
%     [control, state] = catheter.min_contact_(velocity_samples, alpha, state, control, N_x, tip_, disturbances, frictionCoefficient );
% %     [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );
%     P_s = 0; %consider all P_s of all alphas
%     for i = 1: size(alpha,1)
%         w_v = [alpha(i), pi/2 - alpha(i), pi/2];
%         [sigma_mu, f_c, P_s_i] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
%         P_s = P_s + P_s_i;
%     end
%     P_s = P_s / size(alpha,1);
%     
%     k
%     P_s
%     catheter.plot_catheter(state, 'red');
%     pause(0.1);
% end
% 
% %plot final configuration
% catheter.plot_catheter(state, 'g');
% hold on;
% catheter.plot_catheter(state_0, 'blue');
% disp('Enter to proceed');
% pause;

%% evaluate blood flow disturbance
% % % 
alpha_range = [0, 2*pi]';
w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
catheter.velocity_angle_analysis(alpha_range, w_v,  state, control, disturbances, frictionCoefficient);

% velocity_samples = 0.1;
% alpha = 0;
% w_v_alpha = [alpha, pi/2 - alpha, pi/2]';
% externalWrenches = zeros(6,4);
% [F_e] = catheter.compute_external_force(velocity_samples, w_v_alpha, state);
% [f_c, sigma_mu, jacobian] = catheter.contact_force_flow_(state, control, externalWrenches, F_e)

figure(4) % draw P_s w.r.t the angles, given the blood smaple from velocity data
alpha_range = [0, 2*pi]';
w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
catheter.plot_P_s(alpha_range, w_v, state, control, disturbances, frictionCoefficient)