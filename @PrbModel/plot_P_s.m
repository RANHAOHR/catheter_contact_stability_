%%plot the sigma_mu with respect to the blood flow angles
function plot_P_s(obj, alpha_range, w_v, state, control, disturbances, frictionCoefficient)
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
