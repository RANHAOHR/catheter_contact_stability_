function D = material_matrix
    E = 5;
    v = 0.28;

    d_11 = E * (1-v) / ( (1- 2*v) * (1+v) );
    d_12 = E * v / ( (1- 2*v) * (1+v) );

    D_1 = ones(3,3) * d_12 - eye(3,3) * d_12 + eye(3,3) * d_11;

    D_2 = eye(3,3) * (d_11 - d_12) / 2;

    D = zeros(6,6);
    D(1:3,1:3) = D_1;
    D(4:6,4:6) = D_2;

end