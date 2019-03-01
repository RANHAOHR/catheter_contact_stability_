function [F_e] = compute_external_force(obj, velocity_samples, direction_angle, jointAngles)
    
    k_vec = obj.compute_drag_coeffs(direction_angle, jointAngles);

    tip = obj.tip_position(jointAngles);
    positional_difference = zeros(3, obj.numJoints);
    
    g = eye(4);
    for i = 1:1:obj.numJoints
        
        linkConfiguration = obj.linkConfigurations(:, :, i); % the 0 configuration links
        wvec = jointAngles(3*i-2: 3*i,1);
        jointPosition = [0; 0; obj.jointPositions(1, i)];
        g = g * se3rot(wvec, jointPosition, 1);
        linkConfiguration = g * linkConfiguration;
  
        positional_difference(:,i) = linkConfiguration(1:3, 4) - tip;

    end

    Nsample = size(velocity_samples,2);

    F_e = zeros(12, Nsample);
    for j = 1:Nsample
        f_e_all = zeros(3,1);
        tau_e_ = zeros(3,1);
        if velocity_samples(j) < 0
            direction = -direction_angle;
        else
            direction = direction_angle;
        end
        

        externalWrenches = zeros(6, obj.numJoints );
        for i = 1:obj.numJoints
            f_e_ = k_vec(i) * direction * velocity_samples(j) * velocity_samples(j);
            f_e_all = f_e_all + f_e_;
            
            p_it = vector2skewsym(positional_difference(:,i));
            tau_e_ =  tau_e_ + p_it * f_e_; %cross(positional_difference(:,i), f_e_temp);
            externalWrenches(:,i) = [f_e_; tau_e_ ];
        end
        
        
        [F_e(:,j), ~] = obj.external_drag_(jointAngles, externalWrenches);
    end  
    

%end of function    
end