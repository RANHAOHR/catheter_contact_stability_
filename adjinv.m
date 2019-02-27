function Ad_gInv = adjinv(g)
% Find the inverse of Ad_g from g using the equation following eq 2.58 in [1].
    
    R = g(1:3,1:3);
    p = vec_map(g(1:3,4));
    Ad_gInv = [R' -R'*p; zeros(3) R'];
    
end
