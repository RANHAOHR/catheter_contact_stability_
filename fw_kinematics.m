function [exp_cell_, R_cell_, q_init, joint_pts]=fw_kinematics(l_c, l, theta_vec, lineSpec)

    if size(theta_vec, 1) == 12
            jointAngles = theta_vec;
    elseif size(theta_vec, 1) == 3
        jointAngles = reshape(theta_vec,12,1);
    end
    % start Jacobian computation 
    q4 = [0, 0, sum(l(1:4))]';
    q3 = [0, 0, sum(l(1:3))]';
    q2 = [0, 0, sum(l(1:2))]';
    q1 = [0, 0, l(1)]';

    q_init = [q1, q2, q3, q4 ];

    exp_cell_ = cell(1, 4);
    R_cell_ = cell(1,4);
    for i = 1:4
        w_vec = jointAngles(3*i-2:3*i, 1); %check if this is a column vector
        exp_cell_(i) = {compute_exponential_map_(w_vec,q_init(:,i))}; % x_i
        g = cell2mat(exp_cell_(i));
        R_cell_(i) = {g(1:3,1:3)};
    end

    g = eye(4);
    q_ = [q_init, [0;0;l_c]]; %last point is end-effector
    joint_pts = zeros(3,size(q_, 2) ); %include end-effector
    for itr = 1:1:size(q_, 2)
        p = [q_(:,itr);1];
        joint_pt = g * p;
        
        % Draw the catheter
        joint_pts(:,itr) = joint_pt(1:3,1);
        
        if itr < size(q_, 2)
         g = g * cell2mat(exp_cell_(itr));
        end
    end

    plot_catheter(joint_pts, l_c, lineSpec)
    
end