function [k_vec]=compute_drag_coeffs(obj, direction_angle, jointAngles)
    
    outer_diam = obj.outerRadius * 2;
    rho = 1.04 * 10^(-3);
    query = obj.linkLengths(2:end) / outer_diam;
    %approximate according to the book
    dimen_ratio=[1,2,5,10,40,100000];
    cd_ = [0.63,0.68,0.74,0.82,0.98,1.20];
    C_D = interp1(dimen_ratio, cd_, query);

    Aera_0 = outer_diam * obj.linkLengths(2:end); % initial area at zero configuration, linkLengths 5 elements
%     Aera_z = pi * outer_diam^2 /4;
%     Aera_R = zeros(3,5); %for debug
    final_area = zeros(1, obj.numJoints);
    R = eye(3);
    for i = 1:obj.numJoints
        wvec = jointAngles(3*i-2: 3*i,1);
        R = R * so3rot(wvec,1);
        l_i = R * [Aera_0(i),Aera_0(i),0]'; % area in x y z  aidrection after rotation
        l_i = abs(l_i);
%         Aera_R(:,i) = l_i; %the z length
        final_area(i) = compute_area(direction_angle, l_i);
    end

    k_vec = 0.5 * rho *( C_D .* final_area);
end

function [area]=compute_area(direction_angle, area_vec)
    area = sqrt( (area_vec(1) * direction_angle(1))^2 + (area_vec(2) * direction_angle(2))^2 +(area_vec(3) * direction_angle(3))^2 );
end