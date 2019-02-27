function [F_e] = compute_external_force(velocity_samples, direction_angle, exp_cell_,R_cell_, l, g_st)
    k_vec = compute_drag_coeffs(direction_angle, R_cell_, l);
    link_origin_0 = zeros(3,5);
    link_origin_0(3,1) = 0.5 * l(1);
    for i = 2:5
        link_origin_0(3,i) = 0.5 * l(i) + sum(l(1:i-1));
    end
    link_origin = zeros(4,5);
    g = eye(4);
    positional_difference = zeros(3,5);
    for i =1:5
        point_vec = [squeeze(link_origin_0(:,i));1];
        link_origin(:,i) = g * point_vec; % the link origins in spatial frame
        positional_difference(:,i) = link_origin(1:3, i) - g_st(1:3, 4);
        if i < 5
            g = g * cell2mat(exp_cell_(i));
        end
    end

    Nsample = size(velocity_samples,2);

    F_e = zeros(6, Nsample);
    for j = 1:Nsample
        f_e_all = zeros(3,1);
        tau_e_all = zeros(3,1);
        if velocity_samples(j) < 0
            direction = - direction_angle;
        else
            direction = direction_angle;
        end
        for i = 1:5
            f_e_temp = k_vec(i) * direction * velocity_samples(j) * velocity_samples(j);
            f_e_all = f_e_all + f_e_temp;
            
            p_it = vec_map(positional_difference(:,i));
            tau_e_all =  tau_e_all + p_it * f_e_temp; %cross(positional_difference(:,i), f_e_temp);
        end
        F_e(:,j)= [f_e_all; tau_e_all ];

        
    end  
    
    
    
end