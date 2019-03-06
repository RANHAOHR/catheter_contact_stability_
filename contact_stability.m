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

% Calculcate initial configuration
state = [0.5236,0,0,-0.5236,0,0,0.2311,0,0, 0, 0, 0]';
tip_ = catheter.tip_position(state);
control = [0; 0; -1.5; 0; 0; -1.5]; % this can be really far away...

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

[control, state] = catheter.min_contact_( state, control, N_x, tip_, disturbances )
f_c_opt = catheter.contact_force(state, control, disturbances)
sigma_mu_opt = catheter.compute_contact_ratio(f_c_opt)

% catheter.tip_position(state)
% tip_
catheter.plot_catheter(state, 'red');
f_c_0
sigma_mu = catheter.compute_contact_ratio(f_c_0)

disp('Enter to proceed');
pause;

%% Add blood flow disturbance


figure(2)
[velocity_samples] = blood_flow; % get blood flow samples

% velocity_samples = [0:0.1:0.9];
alpha = [0:0.5:2*pi]';
ps_1 = zeros(size(alpha,1),1);
for i = 1: size(alpha,1)
        beta = [alpha(i),pi/2 - alpha(i),pi/2];
        [sigma_mu1, f_c1 ps_1(i)] = catheter.compute_sigma_(velocity_samples, beta, state, control, disturbances, frictionCoefficient );
        i
end
plot(alpha, ps_1, 'b-.','LineWidth',2 );
hold on;

% 
% alpha = [-pi/6:0.01:pi/6]';
% sigma_mu1 = zeros(1,size(alpha,1));
% for j = 1:size(velocity_samples,2)
%     for i = 1: size(alpha,1)
%         beta = [pi/2, alpha(i),pi/2 - alpha(i)];
%         [sigma_mu1(i), f_c1 ps_1] = catheter.compute_sigma_(velocity_samples(j), beta, state, control, disturbances, frictionCoefficient );
%         
%     end
%     plot(alpha, sigma_mu1, 'b-.','LineWidth',2 );
%     hold on;
% 
%     pause;
%     j
% end
% xlabel('Blood flow directional angle $\theta$ /rad','Interpreter','latex');
% ylabel('Contact ratio $\sigma_\mu$ ','Interpreter','latex');
% 

xlabel('Blood flow directional angle $\theta$ /rad','Interpreter','latex');
ylabel('Contact Stability Measure $P_s$ ','Interpreter','latex');


Setfontsize = 15;
set(get(gca,'Xlabel'),'FontSize',Setfontsize);
set(get(gca,'Ylabel'),'FontSize',Setfontsize);
