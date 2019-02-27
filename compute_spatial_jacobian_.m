function [J_sf]=compute_spatial_jacobian_(q_init,jointAngles)
    J_sf = zeros(6,12);
    configurations = zeros(4, 4, 4);
    w_vec = jointAngles(1:3, 1);
    configurations(:, :, 1) = compute_exponential_map_(w_vec, q_init(:,1));    
    
    for i = 1:3
        w_vec = jointAngles(3*(i+1)-2:3*(i+1), 1);
        configurations(:, :, i + 1) = configurations(:, :, i) * compute_exponential_map_(w_vec, q_init(:,(i+1)));
    end

    % Calculate the spatial Jacobian.
    for i = 1 : 4
        w_vec = jointAngles(3*i-2:3*i, 1);
        [wx, wy, wz] = velocity_axes(w_vec);
        
        adjoint = adj(configurations(:, :, i));
        J_sf(:, 3 * i - 2) = adjoint * [-cross(wx,q_init(:,i));wx];
        J_sf(:, 3 * i - 1) = adjoint * [-cross(wy,q_init(:,i));wy];
        J_sf(:, 3 * i) = adjoint * [-cross(wz,q_init(:,i));wz];   
       
    end

end
