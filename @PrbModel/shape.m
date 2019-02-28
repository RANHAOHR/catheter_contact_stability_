function positions = shape(obj, angles)
%
% positions = shape(obj, angles)
%
% Returns locations of all the joints in the spatial frame.
%
% Input:
% angles is a joint angle vector
%
% Output:
% positions is an array of joint positions. Each position is an R^3
% row-vector.
%

% Base, joints, and tip positions
points = [0, obj.jointPositions, obj.length];
% Calculate configurations
configurations = obj.point_configurations(angles, points);
% Get position
positions = configurations(1:3, 4, :);
% Reshape position
positions = reshape(positions, [3, size(configurations, 3)]);

end
