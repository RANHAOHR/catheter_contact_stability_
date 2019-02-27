function [tau_u, B, P]=compute_tau(q_init, jointAngles, u)


   [B, P, J_u, J_b_c1, J_b_c2] = actuation_map(q_init, jointAngles);


   g = [0,0,9.8,0,0,0]';
      
   tau_mg = J_b_c1'* 10^(-4) * g + J_b_c2'* 10^(-4) * g;
      
   tau_u = J_u' * B * P * u + tau_mg;

end