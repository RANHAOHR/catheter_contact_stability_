function L_ = partial_L(x,y,z, phi_)
    
    phi_row1 = diff(phi_, x);
    phi_row2 = diff(phi_, y);
    phi_row3 = diff(phi_, z);
    
    phi_row4 = diff(phi_, z) + diff(phi_, y);
    phi_row5 = diff(phi_, z) + diff(phi_, x);
    phi_row6 = diff(phi_, y) + diff(phi_, z);
    
    L_ = [phi_row1;phi_row2;phi_row3;phi_row4;phi_row5;phi_row6];
  
end