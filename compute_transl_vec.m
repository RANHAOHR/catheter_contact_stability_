function [trans_vec]=compute_transl_vec(q,ew)
 I = eye(3);
 trans_vec = (I - ew)*q;
 
end