function tip_position = tip_position(gst_0, q_init, jointAngles)

    g_st = eye(4);
    for i = 1 : 4
        w_vec = jointAngles(3*i-2:3*i, 1); %check if this is a column vector
        g = compute_exponential_map_(w_vec,q_init(:,i)); % x_i
        g_st = g_st * g; 
    end
    
    g_st = g_st * gst_0;
    tip_position = g_st(1:3,4);

end
