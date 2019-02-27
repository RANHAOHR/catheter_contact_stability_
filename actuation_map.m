function [B,P, J_u,J_b_c1, J_b_c2] = actuation_map(q_init, jointAngles)
    
    gst_c1 = [1,0,0,0;
             0,1,0,0;
             0,0,1,0.044;
             0,0,0,1];
    [J_b_c1,J_u1] = compute_actuator_jacobian_(gst_c1, q_init, jointAngles, 2);%actuator 1 is in the 2rd joint

    gst_c2 = [1,0,0,0;
             0,1,0,0;
             0,0,1,0.0752;
             0,0,0,1];
    [J_b_c2,J_u2]  = compute_actuator_jacobian_(gst_c2, q_init, jointAngles, 3); %actuator 2 is in the 3rd joint

    J_u = [J_u1;J_u2];

    B_0 = [0,3,0]'; %3T homogeneous magnetic field
    B_c = zeros(3,1,2);
    
    B = zeros(3 * 2, 2 * 2);  % Init compact representation of input matrix: actuatorDofs * m, actuatorTorqueDofs * m
    P = zeros(2 * 2, 3 * 2);  % Init mapping from compact to non-compact input    

    nz = 100;
    nx = 30;
    ny = 30; 

    w_vec = jointAngles(1:3, 1); %check if this is a column vector
    g1 = compute_exponential_map_(w_vec,q_init(:,1)); % x_i
    w_vec = jointAngles(4:6, 1); %check if this is a column vector
    g2 = compute_exponential_map_(w_vec,q_init(:,2)); % x_i
    R = g1(1:3,1:3) * g2(1:3,1:3);
        
    for i = 1:2
        turnAreaMatrices = 1e-3 * [nx * 46.6650, ny * 46.6650, nz * 11.251791] .* eye(3);
        B_c(:,:,i) = R * B_0;
        w_vec = jointAngles(7:9, 1); %check if this is a column vector
        g = compute_exponential_map_(w_vec,q_init(:,2)); % x_i       
        R = R * g(1:3,1:3);
        
        % Take out null space of magnetic field via SVD
        try
            [U, S, V] = svd(-vec_map(B_c(:,:,i)) * turnAreaMatrices); % SVD of -B_b * NA (minus sign from reversing order of cross product)
        catch
            warning('Problem with SVD of input matrix');
        end  
        B(3*i-2:3*i, 2*i-1:2*i) = U(1:3, 1:2) * S(1:2, 1:2);   % Calculate the compact representation of B_b
        P(2*i-1:2*i, 3*i-2:3*i) = V(1:3, 1:2)' ;                   
        
    end
end