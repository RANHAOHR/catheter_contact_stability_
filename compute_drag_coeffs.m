function [k_vec]=compute_drag_coeffs(direction_angle, R_cell_, l)
    outer_diam = 3.2 * 10^(-3);
    Aera_0 = outer_diam * l;
%     Aera_z = pi * outer_diam^2 /4;
%     Aera_R = zeros(3,5); %for debug
    final_area = zeros(1,5);
    R = eye(3);
    for i = 1:5
        l_i = R * [Aera_0(i),Aera_0(i),0]';
%         Aera_R(:,i) = l_i; %the z length
        final_area(i) = compute_area(direction_angle, l_i);
        if i < 5
          R = cell2mat(R_cell_(i)) * R;
        end
    end
%     Aera_R
%     final_area

    rho = 1.04 * 10^3;
    query = l / outer_diam;
    %approximate according to the book
    dimen_ratio=[0.5, 1,2,5,10,40,100000];
    cd_ = [0.6,0.63,0.68,0.74,0.82,0.98,1.20];
    C_D = interp1(dimen_ratio, cd_, query);

    k_vec = 0.5 * rho *( C_D .* final_area);
end