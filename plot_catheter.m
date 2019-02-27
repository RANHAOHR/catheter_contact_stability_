function plot_catheter(joint_pts, l_c, lineSpec)

    if nargin == 2
        lineSpec = 'r';
    end
    % Initialize figure
    % figure; 
    hold on; 
    grid on;
    len = l_c;
    axis([-len, len, -len, len, 0, len]);
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    
    % Plot catheter for each joint angle vector
    joint_pts = [[0;0;0], joint_pts];
    plot3(joint_pts(1, :),joint_pts(2, :),joint_pts(3, :), lineSpec, 'LineWidth',2);
    Setfontsize = 15;
    set(get(gca,'Xlabel'),'FontSize',Setfontsize);
    set(get(gca,'Ylabel'),'FontSize',Setfontsize);
    hold off;
end
