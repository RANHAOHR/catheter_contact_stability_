%%plot the sigma_mu with respect to the blood flow angles
function velocity_angle_analysis(obj, alpha_range, w_v, state, control, disturbances, frictionCoefficient)
figure(2)

velocity_samples = [0:0.1:0.9];
alpha = [alpha_range(1):0.1:alpha_range(2)]';
w_v_alpha = w_v(alpha);
sigma_mu1 = zeros(1,size(alpha,1));
for j = 1:size(velocity_samples,2) % draw the sigma w.r.t the angles, given each velocity scale level
    for i = 1: size(alpha,1)
        [sigma_mu1(i), f_c1, ps_1] = obj.compute_sigma_(velocity_samples(j), w_v_alpha(i,:), state, control, disturbances, frictionCoefficient );
        
    end
    x = alpha * 360 / (2*pi);
    plot(x, sigma_mu1, 'b-.','LineWidth',2 );
    hold on;

    pause;
    j
end
ylim([0 0.3]);
xlabel('Blood flow directional angle $\alpha$  /$^\circ$','Interpreter','latex');
ylabel('Contact ratio $\sigma_\mu$ ','Interpreter','latex');

Setfontsize = 15;
set(get(gca,'Xlabel'),'FontSize',Setfontsize);
set(get(gca,'Ylabel'),'FontSize',Setfontsize);

disp('plot y-z plane');
pause;

figure(3) % draw P_s w.r.t the angles, given the blood smaple from veloti data

[velocity_samples] = blood_flow;

alpha = [alpha_range(1):0.1:alpha_range(2)]';
w_v_alpha = w_v(alpha);
ps_1 = zeros(size(alpha,1),1);
for i = 1: size(alpha,1)
        [sigma_mu1, f_c1, ps_1(i)] = obj.compute_sigma_(velocity_samples, w_v_alpha(i,:), state, control, disturbances, frictionCoefficient );
       
end

x = alpha * 360 / (2*pi);
plot(x, ps_1, 'b-.','LineWidth',2 );
ylim([0.5 1]);
hold on;

xlabel('Blood flow directional angle $\alpha$  /$^\circ$','Interpreter','latex');
ylabel('Contact Stability Measure $P_s$ ','Interpreter','latex');


Setfontsize = 15;
set(get(gca,'Xlabel'),'FontSize',Setfontsize);
set(get(gca,'Ylabel'),'FontSize',Setfontsize);
