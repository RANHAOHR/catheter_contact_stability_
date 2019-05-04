function L_ = partial_L(x,y,z, phi_i)
    
    phi_row1 = diff(phi_i, x);
    phi_row2 = diff(phi_i, y);
    phi_row3 = diff(phi_i, z);
    
    phi_row4 = diff(phi_i, z) + diff(phi_i, y);
    phi_row5 = diff(phi_i, z) + diff(phi_i, x);
    phi_row6 = diff(phi_i, y) + diff(phi_i, z);
    
    L_ = [phi_row1;phi_row2;phi_row3;phi_row4;phi_row5;phi_row6];
  
end