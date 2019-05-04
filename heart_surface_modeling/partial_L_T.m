function L_T = partial_L_T(x,y,z, stress)
    
    sigma_xx = stress(1,:);
    sigma_yy = stress(2,:);
    sigma_zz = stress(3,:);
    
    sigma_yz = stress(4,:);
    sigma_xz = stress(5,:);
    sigma_xy = stress(6,:);
    
    phi_row1 = diff(sigma_xx, x) + diff(sigma_xz, z) + diff(sigma_xy, y);
    phi_row2 = diff(sigma_yy, y) + diff(sigma_yz, z) + diff(sigma_xy, x);
    phi_row3 = diff(sigma_zz, z) + diff(sigma_yz, y) + diff(sigma_xz, x);

    L_T = [phi_row1; phi_row2; phi_row3];
  
end