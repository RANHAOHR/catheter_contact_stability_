function [N_theta]=compute_N_theta( theta_vec, l )

    mass_vec = [0.01655,0.2181,0.2893,0.2455,0.2306] * 1.2 * 10^(-3); %kg
    N_theta = zeros(12,1);
    
    for i = 1:3
        [d_theta_x, d_theta_y, d_theta_z] = get_deriv(theta_vec,i);
        N_theta(3*i -2) = ( 0.5 * mass_vec(i+1) + sum(mass_vec(i+2:5)) ) * 9.8 * l(i+1) * d_theta_x;
        N_theta(3*i -1) = ( 0.5 * mass_vec(i+1) + sum(mass_vec(i+2:5)) ) * 9.8 * l(i+1) * d_theta_y;
        N_theta(3*i ) = ( 0.5 * mass_vec(i+1) + sum(mass_vec(i+2:5)) ) * 9.8 * l(i+1) * d_theta_z;
 
    end  
    [d_theta_x, d_theta_y, d_theta_z] = get_deriv(theta_vec, 4);

    N_theta(10) = ( 0.5 * mass_vec(5) ) * 9.8 * l(5) * d_theta_x;
    N_theta(11) = ( 0.5 * mass_vec(5) ) * 9.8 * l(5) * d_theta_y;
    N_theta(12) = ( 0.5 * mass_vec(5) ) * 9.8 * l(5) * d_theta_z;
    
end


function [d_theta_x, d_theta_y, d_theta_z] = get_deriv(theta_vec,i)

    theta_3_2 = theta_vec(1,i)^2 + theta_vec(2,i)^2 +theta_vec(3,i)^2 ;
    eps = 1e-16;
   if theta_3_2 < eps 
        d_theta_x = 0;
        d_theta_y = 0;
        d_theta_z = 0;
   else
    t1 = -sin(sqrt( theta_3_2 ) ) * theta_vec(1,i) * (theta_vec(1,i)^2 + theta_vec(2,i)^2) / sqrt( theta_3_2) +  2*theta_vec(1,i) * cos(sqrt( theta_3_2));
    n1 = t1 * (theta_3_2) - 2 * theta_vec(1,i) * cos(sqrt( theta_3_2)) *( theta_vec(1,i)^2 + theta_vec(2,i)^2);
    n2 = 2 *  theta_vec(1,i) *  theta_vec(3,i)^2;
    d_theta_x = n1 - n2 / theta_3_2;
    
    t2 = -sin(sqrt( theta_3_2 ) ) * theta_vec(2,i) * (theta_vec(1,i)^2 + theta_vec(2,i)^2) / sqrt( theta_3_2 ) +  2*theta_vec(2,i) * cos(sqrt( theta_3_2));
    n3 = t2 * (theta_3_2) - 2 * theta_vec(2,i) * cos(sqrt( theta_3_2)) * (theta_vec(1,i)^2 + theta_vec(2,i)^2);
    n4 = 2 *  theta_vec(2,i) *  theta_vec(3,i)^2;
    d_theta_y = n3 - n4 / theta_3_2;   
  
    t3 = -sin(sqrt( theta_3_2 ) ) * theta_vec(3,i) * (theta_vec(1,i)^2 + theta_vec(2,i)^2) / sqrt(theta_3_2 ) * theta_3_2;
    n5 = t3 - 2  * theta_vec(3,i) * cos(sqrt( theta_3_2)) * (theta_vec(1,i)^2 + theta_vec(2,i)^2);
    n6 = 2 *  theta_vec(3,i) *  (theta_vec(1,i)^2 + theta_vec(2,i)^2);
    d_theta_z = n5 + n6 / theta_3_2;   
   end

end