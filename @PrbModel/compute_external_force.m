function [F_e] = compute_external_force(obj, velocity_samples, direction_angle, jointAngles)
    
    if size(direction_angle,2)>1
        disp('error: use row vector');
        return
    end
    k_vec = obj.compute_drag_coeffs(direction_angle, jointAngles);
    
    jointLinkVecs = zeros(3, obj.numJoints, obj.numJoints);
    linkPositions_ = zeros(3, obj.numJoints);
    jointPositions_ = zeros(3, obj.numJoints);
    
    g = eye(4);
    for i = 1:1:obj.numJoints
        
        linkConfiguration = obj.linkConfigurations(:, :, i); % the 0 configuration links
        wvec = jointAngles(3*i-2: 3*i,1);
        jointPosition = [0; 0; obj.jointPositions(1, i)];
        g = g * se3rot(wvec, jointPosition, 1);
        linkConfiguration = g * linkConfiguration;
        
        linkPositions_(:,i) = linkConfiguration(1:3, 4);
        jointPosition_ = g * [0; 0; obj.jointPositions(1, i);1];    
        jointPositions_(:,i) = jointPosition_(1:3);

    end

    for j = 1:obj.numJoints
      jointPosition_ = repmat(jointPositions_(:,j), [1,obj.numJoints]);
      jointLinkVecs(:,:,j) = linkPositions_ - jointPosition_; %should be a 3 by 4 by 4, each (:,:,j) represent the link to the joint positional differences 
    end
    
    Nsample = size(velocity_samples,2);
    
    direction = repmat(direction_angle, [1,Nsample]);
    [~, ind] = find(velocity_samples < 0); 
    direction(ind) = -direction(ind); % flip the directions
    v_ = velocity_samples .^2;
    v_ = repmat(v_, [3,1]); % for elementwise multiplication
    v2_vec =  direction .* v_;
    
    f_e_ = zeros(3, Nsample, obj.numJoints );
    
    externalWrenches = zeros(6, Nsample, obj.numJoints );
    for i = 1:obj.numJoints
        f_e_(:,:,i) = k_vec(i) .* v2_vec; % for all forces, this may change if the direction in each link is different
    end

    for i = 1:obj.numJoints
        
        tau_e_ = zeros(3, Nsample);
        link_diffs = jointLinkVecs(:,:,i);       
        for j = 1:obj.numJoints
            p_it = vector2skewsym(link_diffs(:,j));
            tau_e_ = p_it * f_e_(:,:,j) + tau_e_; %the torque for one link of all blood samples
        end

        externalWrenches(:,:,i) = [f_e_(:,:,i); tau_e_];
    end
    
    [F_e, ~] = obj.external_drag_(jointAngles, externalWrenches);
 
%end of function    
end