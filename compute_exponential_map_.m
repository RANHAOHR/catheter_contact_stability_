function [exp_map]=compute_exponential_map_(w_vec,q)
    exp_map = zeros(4);
    exp_map(4,4) = 1;
    
    expoential_w = compute_ew_(w_vec); % indice 0,1,2
    vec = compute_transl_vec(q,expoential_w);
    
    for i = 1:3
        temp = squeeze(expoential_w(i,:));
        exp_map(i,1:3) = temp;
        exp_map(i,4) = vec(i);
    end
end