function x = tip_position(obj, q)
%
% x = tipPosition(q)
%
% Calculate tip position, x, from given joint angles, q.
%
% Input:
% q is the joint angles.
%
% Output:
% x is the catheter's tip position.
%

% Calculate tip configurations
configurations = obj.point_configurations(q, [obj.length]);
% Get position
x = configurations(1:3, 4, 1);

end
