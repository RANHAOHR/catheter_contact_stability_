function Ad_g = adj(g)
% Find Ad_g from g using eq 2.58 in [1].

    R = g(1:3,1:3);
    p = vec_map(g(1:3,4));
    Ad_g = [R p*R; zeros(3) R];

end
