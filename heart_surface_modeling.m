% Clear workspace.
clear;
close all;
clc;
set(0,'defaultfigurecolor','w');

%test in 1D
% min = 1;
% max = 5;
% 
% Nsmaple = 50;
% x = min + rand(1,Nsmaple)*(max - min);
% % x= sort(x);
% y = min + rand(1,Nsmaple)*(max - min);
% 
% z = sin(x+0.5*y);% + rand(1,20);
% 
% samples = [x;y];
% 
% figure(1);
% scatter3(x, y, z, 'r.');
% hold on;
% query_ = [1.0,3.5];
% p_x = [1, query_(1), query_(2), query_(1) * query_(2), query_(1)^2, query_(2)^2];
% kernel_size = 6;
% P_I = zeros(Nsmaple, kernel_size);
% W_I = zeros(Nsmaple);
% for i = 1:Nsmaple
%     P_I(i,:) = [1, x(i), y(i), x(i) * y(i), x(i)^2, y(i)^2];
%     d = norm(query_ - samples(:,i),2);
%     W_I(i,i) = exp(-5 * d^2); % should be good enough to exagerate differences, coefficient should not be to big
% 
% end
% W_I = W_I ./ norm(W_I,1);
% P_I = P_I';
% A = P_I * W_I * P_I';
% 
% B = P_I * W_I;
% a_x = A \ B * z';
% 
% query_y = p_x * a_x;
% 
% scatter3(query_(1), query_(2), query_y, 'b.');

%surface modeling
N = 20;
phi_ = kernel_3d(N)


xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');