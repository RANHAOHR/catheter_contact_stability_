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
% state = [-0.624852668813324;-0.0542563518171362;-0.0525068370193255;-0.676749434104250;-0.0587686506749256;-0.0524679635706560;0.313943952439617;-0.336741809028734;0.0498084449473635;-0.0130426648031853;-0.0483865926995447;2.42792650671112e-08];
% tip_ = catheter.tip_position(state);
% control = [-0.0813856046291852;-0.213973639206613;0.273749605707889;0.112216434314529;0.157425830943400;0.0853017218434574]; % this can be really far away...

% 
% state = [-0.0737020339466517;-0.0218040664877389;0.557950048047828;-0.255235186705955;0.337803859876394;0.557141793151812;-0.179753774914578;-0.0895241380738602;1.55912523055963;-0.240881919228409;-0.175779886865597;1.43491696547212e-05];
% tip_ = catheter.tip_position(state);
% control = [-0.318923394305813;0.247555477222941;1.35576833186285;0.0412604294881832;-1.31774443872355;1.48556110031736];

state = [0.0725, -0.0514, 1.1160, 0.1259, -0.2842, 1.1149, -0.2702, -0.1743, -1.5707, -0.2305, 0.2968, 0.0]';
tip_ = catheter.tip_position(state);
control = [0.6544, -2.0530, -1.5000, -0.08723, 1.5228, -1.9458]'; % this can be really far away...

origin = [0; 0; tip_(3)] %% just to gurantee on the surface plane now...
orientation = so3rot([1; 0; 0], pi);
frictionCoefficient = 0.2;
plane = Plane(origin, orientation, frictionCoefficient);
catheter.set_surface(plane);
tip_ = catheter.tip_position(state);

% Print results
figure(1);
disp('Original contact force without blood flow = ');
f_c_0 = catheter.contact_force(state, control, disturbances)
% catheter.tip_position(state0);
catheter.plot_catheter(state, 'blue');
pause;

%% evaluate blood flow disturbance

% alpha_range = [0, 2*pi]';
% w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
% catheter.velocity_angle_analysis(alpha_range, w_v,  state, control, disturbances, frictionCoefficient);

[velocity_samples] = blood_flow;

alpha = pi;
w_v = [alpha, pi/2 - alpha, pi/2];

[~, ~, P_s] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
display('Initial P_s under the given blood flow angle: ');
P_s
pause;
%% start contact force control given blood flow

velocity_samples = 0.9;

direction_angle = [cos(w_v(1)),cos(w_v(2)),cos(w_v(3))]';

[F_e] = catheter.compute_external_force(velocity_samples, direction_angle, state); %get the external motion caused by the blood flow
Nsample = size(F_e, 2);

[state, exitflag, lambdas, hessian] = catheter.min_potential_energy_conf_const( state, control, disturbances, [], tip_, options);

[sigma_mu, f_c_0, P_s] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
display('Contact force under blood flow with 0.9m/s = ');
f_c_0

k_f = 0.0001;
k_c = 0.0002;
dfc_x =  - k_f * f_c_0(1) / sqrt(f_c_0(1)^2 + f_c_0(2)^2);
dfc_y =  - k_f * f_c_0(2) /  sqrt(f_c_0(1)^2 + f_c_0(2)^2);
if f_c_0(3) < 0.0
    dfc_z = k_c;
else
    dfc_z = 0.0;
end

dfc = [dfc_x, dfc_y, dfc_z]';

disp('The initial tip position = ');
tip_ = catheter.tip_position(state)
disp('The initial joint angle = ');
state
pause;

for i = 1:500
    
[J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(state, control, disturbances); %endEffectorNullspace here is the transpose fo the null space...

[J_cu, J_ctheta, J_cq] = catheter.compute_contact_jacbobian(state, control, F_e, disturbances);

N_k = eye(12) - pinv(J_e) * J_e;

J_test = J_ctheta * N_k * J_cq + J_cu;
du = pinv(J_test) * dfc;

control = control + du;
state = state + N_k * J_cq * du;
[f_c, ~, jacobian] = catheter.contact_force_flow_(state, control, disturbances, F_e);
f_c
tip_ = catheter.tip_position(state)
catheter.plot_catheter(state, 'red');
i
pause;
end

display('Optimized contact force under blood flow with 0.9m/s: ');
f_c
display('The new state: ');
state

pause;
%% test
[velocity_samples] = blood_flow;
[sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
display('Optimized P_s under the given blood flow angle: ');
P_s
% figure(5);
% alpha_range = [0, 2*pi]';
% w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
% catheter.velocity_angle_analysis(alpha_range, w_v,  state, control, disturbances, frictionCoefficient);