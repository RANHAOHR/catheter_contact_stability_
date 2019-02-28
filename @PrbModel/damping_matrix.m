function C = damping_matrix(obj, q)
%
% C = damping_matrix(q)
%
% This function returns the damping matrix.
%
% Input:
% q is the joint angle vector.
%
% Output:
% C is the fluid damping matrix.
%

remainingLengths = [zeros(obj.numJoints, 1), triu(ones(obj.numJoints, obj.numJoints), 0)] * obj.linkLengths';
C = obj.damping * diag(kron(remainingLengths, ones(3, 1))) / obj.length;

end
