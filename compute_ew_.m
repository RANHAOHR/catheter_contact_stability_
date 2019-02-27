function [ew]=compute_ew_(w_vec)
    eps = 1e-16;
   if norm(w_vec, 2) < eps %avoid singularity
        ew = eye(3);
        % Otherwise calculate rotation matrix
   else
        W_ = vec_map(w_vec);
        ew = eye(3) + W_/sqrt(w_vec' * w_vec) * sin( sqrt(w_vec' * w_vec) ) + W_*W_/sqrt(w_vec' * w_vec)^2 * ( 1-cos(sqrt(w_vec' * w_vec)) );
    end
end