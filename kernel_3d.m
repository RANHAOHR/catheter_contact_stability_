function phi_ = kernel_3d(x, y, z, N)

    p_x = [1];

    kernel = [x, y, z];
    i = 1;
    while i < N
        p_1 = x * kernel;
        p_2 = y * kernel;
        p_3 = z * kernel;

        p_ = [p_1, p_2, p_3];

        kernel = unique(p_);
        p_x = [p_x, kernel];

        if size(p_x, 2) > N
            break;
        end
        i = i+1;

    end
    phi_ = p_x(1:N);
    
    L_u = partial_L(x,y,z, phi_);
  
end