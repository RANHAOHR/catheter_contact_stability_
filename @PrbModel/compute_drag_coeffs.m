function [k_vec]=compute_drag_coeffs(obj, direction_angle, jointAngles)
    
    outer_diam = obj.outerRadius * 2;
    rho = 1.04 * 10^(-3);
    query = obj.linkLengths(2:end) / outer_diam;
    %approximate according to the book
    dimen_ratio=[1,2,5,10,40,100000];
    cd_ = [0.63,0.68,0.74,0.82,0.98,1.20];
    C_D = interp1(dimen_ratio, cd_, query);

    Aera_0 = outer_diam * obj.linkLengths(2:end); % initial area at zero configuration, linkLengths 5 elements

    final_area = zeros(1, obj.numJoints);
    R = eye(3);
    for i = 1:obj.numJoints
        wvec = jointAngles(3*i-2: 3*i,1);
        R = R * so3rot(wvec,1);
        l_i = R * [0,0,Aera_0(i)]'; % area in x y z  aidrection after rotation

        l_angle = l_i / norm(l_i);
        cos_theta = l_angle' * direction_angle;
        final_area(i) = Aera_0(i) * (1 - cos_theta^2);
    end

%     final_area %debug
    k_vec = 0.5 * rho *( C_D .* final_area);
end