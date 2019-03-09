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
% state0 = [0.5236,0,0,-0.5236,0,0,0.2311,0,0, 0, 0, 0]';
% tip_ = catheter.tip_position(state0);
% control0 = [0; 0; 3; 0; 0; -3]; % this can be really far away...

state0 = [0.0725, -0.0514, 1.1160, 0.1259, -0.2842, 1.1149, -0.2702, -0.1743, -1.5707, -0.2305, 0.2968, 0.0]';
tip_ = catheter.tip_position(state0);
control0 = [0.6544, -2.0530, -1.5000, -0.08723, 1.5228, -1.9458]'; % this can be really far away...

% state0 = [-0.5236,0,0,0.5236,0,0,-0.2311,0,0, 0, 0, 0]';
% tip_ = catheter.tip_position(state0);
% control0 = [0.0; 0; 1; 0.0; -0.1; 1.5]; % this can be really far away...

% state0 = [-0.0374;0.04693;0.9673;-0.05602;0.3445;0.96563;0.02618;0.08565;0.6394;-0.3222;-0.2276;-0.00058];
% tip_ = catheter.tip_position(state0);
% control0 = [0.1427;-0.1476;1.0;0.0292;-0.3845;1.4167]; % this can be really far away...

% state0 = [-0.701265863284934;-0.0426838953283018;0.0880352838146160;-0.746305331659681;-0.0268660573931596;0.0879736506997280;0.385266497986597;-0.178974452627024;-0.0899008546983464;-0.0191629314890400;-0.0348938117840623;-6.99294332058007e-09];
% tip_ = catheter.tip_position(state0);
% control0 = [-0.0134330693468157;-0.280794893099882;0.257841708964121;0.0636627849085847;0.178443886519327;0.0920929306759346]; % this can be really far away...

% state0 = [-0.0737020339466517;-0.0218040664877389;0.557950048047828;-0.255235186705955;0.337803859876394;0.557141793151812;-0.179753774914578;-0.0895241380738602;1.55912523055963;-0.240881919228409;-0.175779886865597;1.43491696547212e-05];
% tip_ = catheter.tip_position(state0);
% control0 = [-0.318923394305813;0.247555477222941;1.35576833186285;0.0412604294881832;-1.31774443872355;1.48556110031736]; % this can be really far away...

% Create a plane
origin = [0; 0; tip_(3)]; %% just to gurantee on the surface plane now...
orientation = so3rot([1; 0; 0], pi);
frictionCoefficient = 0.2;

plane = Plane(origin, orientation, frictionCoefficient);
catheter.set_surface(plane);

% Print results
figure(1);
disp('contact_force = ');
f_c_0 = catheter.contact_force(state0, control0, disturbances)
sigma_mu_0 = catheter.compute_contact_ratio(f_c_0)
catheter.tip_position(state0)
catheter.plot_catheter(state0, 'blue');
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

state = state0;
control = control0;
%% start contact force control given blood flow
figure(2);
[velocity_samples] = blood_flow;

alpha = pi;
% w_v = [alpha, pi/2 - alpha, pi/2];
% 
% [sigma_mu, f_c, P_s] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
% pause;


% alpha = [3.14:0.1:4.1733]';
disp('Initial P_s:');
P_s = 0; %consider all P_s of all alphas
for i = 1: size(alpha,1)
    w_v = [alpha(i), pi/2 - alpha(i), pi/2];
    [sigma_mu, f_c, P_s_k] = catheter.compute_sigma_(velocity_samples, w_v, state, control, disturbances, frictionCoefficient );
    P_s = P_s + P_s_k;
end
P_s = P_s / size(alpha,1)

alpha_range = [0, 2*pi]';
w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
catheter.plot_P_s(alpha_range, w_v, state, control, disturbances, frictionCoefficient);
pause;

figure(3);

Ntrial = 4;

states = zeros(3 * catheter.get_num_joints(), Ntrial+1);
states(:, 1) = state;
controls = zeros(3 * catheter.get_num_coilsets(), Ntrial+1);
controls(:, 1) = control;

for k = 1:Ntrial
    tip_ = catheter.tip_position(state);
    [J_e, endEffectorNullspace, J_q, surfaceJacobian] = catheter.jacobian(states(:,k), controls(:,k), disturbances); %endEffectorNullspace here is the transpose fo the null space...

    J_x = J_e * J_q;
    N_x = null(J_x);

    [controls(:,k+1), states(:,k+1)] = catheter.min_contact_(velocity_samples, alpha, states(:,k), controls(:,k), N_x, tip_, disturbances, frictionCoefficient );

    P_s = 0; %consider all P_s of all alphas
    for i = 1: size(alpha,1)
        w_v = [alpha(i), pi/2 - alpha(i), pi/2];
        [~, ~, P_s_k] = catheter.compute_sigma_(velocity_samples, w_v, states(:,k+1), controls(:,k+1), disturbances, frictionCoefficient );
        P_s = P_s + P_s_k;
    end
    P_s = P_s / size(alpha,1);
    
    k
    P_s
    
%     f_c = catheter.contact_force(states(:,k+1), controls(:,k+1), disturbances)
%     sigma_mu = catheter.compute_contact_ratio(f_c)
    pause;
end

figure(4)
%plot final configuration
catheter.plot_catheter(states(:,1), 'blue');
hold on;
for k = 2:Ntrial-1
    catheter.plot_catheter(states(:,k), 'red');
    hold on;
end
catheter.plot_catheter(states(:,end), 'g');

f_c = catheter.contact_force(states(:,end), controls(:,end), disturbances)
sigma_mu = catheter.compute_contact_ratio(f_c)

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

figure(5) % draw P_s w.r.t the angles, given the blood smaple from velocity data
alpha_range = [0, 2*pi]';
w_v = @(alpha)[alpha, pi/2 - alpha, pi/2*ones(size(alpha,1), 1)];
catheter.plot_P_s(alpha_range, w_v, states(:,1), controls(:,1), disturbances, frictionCoefficient);
hold on;
for k = 2:Ntrial-1
    catheter.plot_P_s(alpha_range, w_v, states(:,k), controls(:,k), disturbances, frictionCoefficient, 'r-.');
    hold on;
end
catheter.plot_P_s(alpha_range, w_v, states(:,end), controls(:,end), disturbances, frictionCoefficient, 'g-.');