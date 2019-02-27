clear all;
close all;
clc;
set(0,'defaultfigurecolor','w');

% %compute targte configuration
l_c = 0.104; %m
l_0 = 0.0996;
l = [0.01655,0.2181,0.2893,0.2455,0.2306] * l_c;
% joint_pts = zeros(2,5); %(y,z)
% theta_vec = zeros(1,4);
% joint_pts(:,1) = [0,l(1)]';
% theta_vec(1) = pi/6;
% joint_pts(:,2) = joint_pts(:,1) + l(2) *[-sin(theta_vec(1)), cos(theta_vec(1))]';
% theta_vec(2)= 0;
% joint_pts(:,3) = joint_pts(:,2) + l(3) *[-sin(theta_vec(2)), cos(theta_vec(2))]';
% theta_vec(3) = -abs(asin(joint_pts(1,3) / (l(4) + l(5))));
% joint_pts(:,4) = joint_pts(:,3) + l(4) *[-sin(theta_vec(3)), cos(theta_vec(3))]';
% theta_vec(4) = theta_vec(3); %didn't times rotation here so..
% joint_pts(:,5) = joint_pts(:,4) + l(5) *[-sin(theta_vec(4)), cos(theta_vec(4))]'; %end effector
% 
% figure(4)
% plot([0 joint_pts(1,1)],[0 joint_pts(2,1)],'k','LineWidth',2);
% for i = 1:4
%     hold on;
%     plot([joint_pts(1,i) joint_pts(1,i+1)],[joint_pts(2,i) joint_pts(2,i+1)],'k','LineWidth',2);
% end
% axis([-0.1 0 0 0.1]);

%0.157079632679490 0.071723643759686
theta_vec1 = [[0.5236;0;0], [-0.5236;0;0], [-0.2311;0;0], [0;0;0]];
theta_vec2 = [[-0.5236;0;0], [0.5236;0;0], [0.2311;0;0], [0;0;0]];
theta_vec3 = [[0.5236;0;3.14], [0;0.2211;0], [0;0.2211;0], [0;0.2111;0]];
theta_vec4 = [[0.5236;0;-3.14], [0;-0.2211;0], [0;-0.2211;0], [0;-0.2111;0]];
theta_vec5 = [[0.5236;0;1.57], [-0.2211;0.2211;0], [-0.2211;0.2211;0], [-0.2211;0.2211;0]];
% theta_vec6 = [[0.5236;0;-1.57], [-0.2211;-0.2211;0], [-0.2211;-0.2211;0], [-0.2211;-0.2211;0]];
theta_vec7 = [[0.5236;-1.57;-4.71], [0;-0.2311;0], [-0.2311;-0.2311;0], [0;-0.2311;0]];
figure(4)
% [exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec1, 'r');
% [exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec2, 'r');
[exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec3, 'k');
% [exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec4, 'k');
[exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec5, 'b');
% % [exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec6, 'b');
[exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec7, 'g');
% [exp_cell_, R_cell_, q_init]=fw_kinematics(l_c, l, theta_vec8, 'g');

beta = pi/2;
direction_angle = [cos(beta),cos(pi/2 - beta),0]';
% direction_angle = [0,1,0]';

% velocity_samples = [0:0.01:0.9];
% velocity_samples = 0;
[velocity_samples] = blood_flow(2);

[sigma_mu1, ps1,  P_s1, P_s2, P_s3, P_s4]=catheter_kinematics(direction_angle, velocity_samples, theta_vec1,[0;-0.2;1.06]);
% [sigma_mu2, ps2, P_s21, P_s22, P_s23, P_s24]=catheter_kinematics(direction_angle, velocity_samples, theta_vec2, [0;0.2;1.06]);
% [sigma_mu3, ps3, P_s31, P_s32, P_s33, P_s34]=catheter_kinematics(direction_angle, velocity_samples, theta_vec3, [0.19;0.01;1.06]);
% [sigma_mu5, ps5, P_s51, P_s52, P_s53, P_s54]=catheter_kinematics(direction_angle, velocity_samples, theta_vec5, [0.16;-0.01;1.0]);
% [sigma_mu7, ps7, P_s71, P_s72, P_s73, P_s74]=catheter_kinematics(direction_angle, velocity_samples, theta_vec7, [-0.19;0.03;1.06]); %[-0.19;0.01;1.0]

x = [1:1:10000];
% plot(x, sigma_mu1, 'Color', [0, 0.4470, 0.7410], '.','LineWidth',2 );
% hold on;
% plot(x, sigma_mu2, 'k.','LineWidth',2 );
% hold on;
% plot(x, sigma_mu3, 'r.','LineWidth',2 );
% hold on;
% plot(x, sigma_mu4, 'g.','LineWidth',2 );

% plot(x, sigma_mu1, x,sigma_mu2, x, sigma_mu3, x, sigma_mu5, x,sigma_mu7, x,sigma_mu8,  'LineWidth',2 );
% legend('Configuration 1','Configuration 2','Configuration 3', 'Configuration 5','Configuration 7','Configuration 8' );
% % figure(6) % draw the curve of sigma_relative to value of velocity
% plot(velocity_samples, sigma_mu1, 'b-.','LineWidth',2 );
% hold on;
% plot(velocity_samples, sigma_mu2, 'k-.','LineWidth',2 );
% 
% % [~,ind] = sigma_mu1(0);%  min(sigma_mu1);
% % min_x = velocity_samples(ind);
% 
% xlabel('Blood flow velocity $||v||$ m/s','Interpreter','latex');
% ylabel('Contact ratio $\sigma_\mu$ ','Interpreter','latex');
% legend('Configuration 1','Configuration 2');
% Setfontsize = 15;
% set(get(gca,'Xlabel'),'FontSize',Setfontsize);
% set(get(gca,'Ylabel'),'FontSize',Setfontsize);



% disp('proceed...');
% pause;
% velocity_samples = [0:0.1:0.9];
% % velocity_samples = 0.8;
% % [velocity_samples] = blood_flow(2);
% theta = [0:0.01:pi]';
% % theta = [-pi/6:0.01:pi/6]';
% Ntheta = size(theta,1);
% direction_angle = [cos(theta),sin(theta),zeros(Ntheta,1)]';
% % direction_angle = [zeros(Ntheta,1),cos(theta),sin(theta)]';
% 
% % direction_angle = [0,1,0]';
% sigma_mu1 = zeros(1,Ntheta);
% sigma_mu2 = zeros(1,Ntheta);
% figure(6)
% for j = 1:size(velocity_samples,2)
%     for i = 1:Ntheta
%         direction = squeeze(direction_angle(:,i));
%         [sigma_mu1(i), P_s11, P_s12, P_s13, ps1]=catheter_kinematics(direction, velocity_samples(j), theta_vec1);
%         [sigma_mu2(i), P_s21, P_s22, P_s23, ps2]=catheter_kinematics(direction, velocity_samples(j), theta_vec2);
%     end
% 
%     plot(theta, sigma_mu1, 'b-.','LineWidth',2 );
%     hold on;
%     plot(theta, sigma_mu2, 'k-.','LineWidth',2 );
%     hold on;
%     pause;
%     j
% end
% 
% 
% xlabel('Blood flow directional angle $\theta$ /rad','Interpreter','latex');
% ylabel('Contact ratio $\sigma_\mu$ ','Interpreter','latex');
% 
% 
% legend('Configuration 1','Configuration 2');
% Setfontsize = 15;
% set(get(gca,'Xlabel'),'FontSize',Setfontsize);
% set(get(gca,'Ylabel'),'FontSize',Setfontsize);
% 
