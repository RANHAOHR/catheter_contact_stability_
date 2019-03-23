%% Initialize

% Clear workspace.
clear;
close all;
clc;
set(0,'defaultfigurecolor','w');
% Initialize PrbModel
% surface = plane;
% catheter.set_surface(plane); % set the plane with respect to the tip

% Create a plane
origin = [0; 0; 139];
orientation = so3rot([1; 0; 0], -0.9);
frictionCoefficient = 0.1;
plane = Plane(origin, orientation, frictionCoefficient);

% Initialize PrbModel
surface = plane;
numJoints = 4;
parameterFile = './data/parameters_2coilsets.yml';
catheter = PrbModel(parameterFile, numJoints);
catheter.set_surface(plane);

%% Find quasi-static configuration


% Calculcate initial configuration
state0 = 0.1 * randn(3 * catheter.get_num_joints(), 1);
control0 = [0.05; 0; 0.15; 0.05; 0; 0.15];
disturbances = zeros(6, catheter.get_num_joints());
options = optimoptions('fmincon', 'MaxFunEvals', 1000);
[state, hessian, lambda, exitflag] = catheter.min_potential_energy_conf(...
    state0, control0, disturbances, [], options);

% Print results
disp('contact_force = ');
disp(catheter.contact_force(state, control0, disturbances));
disp('tip_position = ');
disp(catheter.tip_position(state));
disp('distance = ');
disp(surface.distance(catheter.tip_position(state)));

% Calculcate initial configurations
% 
state = [-0.624852668813324;-0.0542563518171362;-0.0525068370193255;-0.676749434104250;-0.0587686506749256;-0.0524679635706560;0.313943952439617;-0.336741809028734;0.0498084449473635;-0.0130426648031853;-0.0483865926995447;2.42792650671112e-08];
tip_ = catheter.tip_position(state);
control0 = [-0.0813856046291852;-0.213973639206613;0.273749605707889;0.112216434314529;0.157425830943400;0.0853017218434574]; % this can be really far away...

% state = [0.0725, -0.0514, 1.1160, 0.1259, -0.2842, 1.1149, -0.2702, -0.1743, -1.5707, -0.2305, 0.2968, 0.0]';
% tip_ = catheter.tip_position(state);
% control0 = [0.6544, -2.0530, -1.5000, -0.08723, 1.5228, -1.9458]'; % this can be really far away...
% 
% origin = [0; 0; tip_(3)] %% just to gurantee on the surface plane now...
% orientation = so3rot([1; 0; 0], pi);
frictionCoefficient = 0.2;
plane = Plane(origin, orientation, frictionCoefficient);
catheter.set_surface(plane);

tip_ = catheter.tip_position(state)
% Print results
figure(1);
disp('contact_force = ');
f_c_0 = catheter.contact_force(state, control0, disturbances)
sigma_mu_0 = catheter.compute_contact_ratio(f_c_0)
% catheter.tip_position(state0);
catheter.plot_catheter(state, 'blue');
pause;

%% start contact force searching

% [J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state0, control0, disturbances); %endEffectorNullspace here is the transpose fo the null space...
% J_x = J_e * J_q;
% N_x = null(J_x);
% 
% velocity_samples = 0;
% alpha = 0;
% 
% [control, state] = catheter.min_contact_(velocity_samples, alpha, state0, control0, N_x, tip_, disturbances, frictionCoefficient );
% state
% control
% [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, alpha, state, control, disturbances, frictionCoefficient );
% catheter.plot_catheter(state, 'red');

control = control0;
%% start contact force control given blood flow
figure(2);
% [velocity_samples] = blood_flow;
velocity_samples = 0;

alpha = 0;
w_v = [alpha, pi/2 - alpha, pi/2];

[sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
sigma_mu
f_c

direction_angle = [cos(w_v(1)),cos(w_v(2)),cos(w_v(3))]';

[F_e] = catheter.compute_external_force(velocity_samples, direction_angle, state); %get the external motion caused by the blood flow
Nsample = size(F_e, 2);

[state, exitflag, lambdas, hessian] = catheter.min_potential_energy_conf_const( state, control, disturbances, [], tip_, options)

dfc1 = 0.0001* f_c(1) / sqrt(f_c(1)^2 + f_c(2)^2);
dfc2 = 0.0001 * f_c(2) / sqrt(f_c(1)^2 + f_c(2)^2);
dfc = [dfc1, dfc2 , 0.0]'
for i = 1:200
    
tip_ = catheter.tip_position(state);
[J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances); %endEffectorNullspace here is the transpose fo the null space...

[f_c, ~, ~] = catheter.contact_force_flow_(state, control, disturbances, F_e);
[J_cu, J_ctheta, J_cq] = catheter.compute_contact_jacbobian(state, control, F_e, J_e, f_c, disturbances);

N_k = eye(12) - pinv(J_e) * J_e;
% dv = zeros(1,size(N_k, 2))';
% dv(1) = 0.001;

dtheta = N_k;

tip_ = catheter.tip_position(state)

J_cq
J_q
J_test = J_ctheta * dtheta * J_cq + J_cu;
du = pinv(J_test) * dfc 

% tip_ = catheter.tip_position(state + J_cq * du)
% dtheta = J_cq * du;
control = control + du;
state = state + dtheta * J_cq * du
[f_c, sigma_mu, jacobian] = catheter.contact_force_flow_(state, control, disturbances, F_e);
f_c
tip_ = catheter.tip_position(state)
catheter.plot_catheter(state, 'red');
pause;
end

state
[state, exitflag, lambdas, hessian] = catheter.min_potential_energy_conf_const( state, control, disturbances, [], tip_, options);
