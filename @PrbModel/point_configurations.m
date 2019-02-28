function configurations = point_configurations(obj, angles, points)
%
% configurations = point_configurations(angles, points)
%
% This function returns configurations of the points along the body
% of the catheter given joint angles.
%
% Input:
% angles is a deflection angle vector.
% points is an array of points along the length of the catheter. The points
% are assumed to be in an increasing order.
%
% Output:
% configurations is an array of configrurations of the points given
% angles. Each element of the array is an elements of SE(3).
%

% Sort points in ascending order.
[sortedPoints, indices] = sort(points);

% Initialize configurations
numPoints = size(sortedPoints, 2);
sortedConfigurations = zeros(4, 4, numPoints);
configuration = eye(4);
% Get link numbers the points are on.
[links, ~] = obj.get_links(sortedPoints);
joints = links;
lastJoint = joints(1, end);
% Forward kinematics initialization
w0 = obj.jointDirections;
p0 = obj.jointPositions;
g = eye(4);
% point index
iPoint = 1;

% Points on the first link do not move
while joints(1, iPoint) == 0
    % Store point configuration
    configuration(3, 4) = sortedPoints(1, iPoint);
    sortedConfigurations(:, :, iPoint) = configuration;
    % Update counter
    iPoint = iPoint + 1;
end

% Calculate configurations of the rest of the points
for i = 1 : 1 : lastJoint
    % Forward kinematics
    w = angles(3 * i - 2 : 3 * i, 1);
    p = [0; 0; p0(1, i)];
    g = g * se3rot(w, p, 1);
    
    % If the point is on this link, calculate and store its configuration
    while iPoint <= numPoints && joints(1, iPoint) == i
        % Calculate and store point configuration
        configuration(3, 4) = sortedPoints(1, iPoint);
        sortedConfigurations(:, :, iPoint) = g * configuration;
        % Update counter
        iPoint = iPoint + 1;
    end
end

% Unsort the configuration array back to the same order as the inputs.
configurations(:, :, indices) = sortedConfigurations;

end
