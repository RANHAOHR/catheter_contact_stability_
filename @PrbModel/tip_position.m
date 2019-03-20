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

% Initialize variables.
configurations = zeros(4, 4, obj.numJoints + 1);
configurations(:, :, 1) = eye(4);

% Calculate the configuration of the joints.
for i = 1 : 1 : obj.numJoints
    configurations(:, :, i + 1) = configurations(:, :, i) * ...
        se3rot(q(3*i-2:3*i, 1), [0; 0; obj.jointPositions(1, i)], 1); 
end

% Calculate tip position at the current configuration.
x = configurations(:, :, end) * [0; 0; obj.length; 1];
x = x(1:3);
% % Calculate tip configurations
% configurations = obj.point_configurations(q, [obj.length]);
% % Get position
% x = configurations(1:3, 4, 1);

end
