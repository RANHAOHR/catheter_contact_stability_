function [linkNumbers, distances] = get_links(obj, points)
%
% [link_numbers, distances] = get_links(obj, points)
%
% This function returns a vector that contains link number which
% points along the length of the catheter is on.
%
% Input:
% points is a list of points along the body of the catheter.
%
% Outputs:
% links is a list of link numbers the points are on.
% distances is a list of distances between the points and the
%  center of masses of the links the points are on.
%
    
% Initialize stuff
linkNumbers = zeros(1, length(points));
distances = zeros(1, length(points));
numPoints = size(points, 2);
sumLengths = tril(ones(obj.numJoints+1)) * obj.linkLengths';

% Loop over points (Lopping the other way is better, I'll do that
% later).
for i = 1:1:numPoints

    % If the point is on the first link.
    if points(1, i) <= obj.jointPositions(1, 1)
        % Store the link the point is on
        linkNumbers(1, i) = 0;
    
    % If the point is beyond the tip, it is considered as being on the last link.
    elseif points(i) > obj.jointPositions(1, end)
        % Store the link the point is on
        linkNumbers(1, i) = obj.numJoints;
        % Calculate distance from point to link center of mass
        distances(1, i) = points(1, i) - obj.jointPositions(1, end);
    end

    for j = 1:1:obj.numJoints-1

        % If point is on link j+1
        if points(1, i) > obj.jointPositions(1, j) && points(1, i) <= obj.jointPositions(1, j+1)
            % Store the link the point is on
            linkNumbers(1, i) = j;
            % Calculate distance from point to link center of mass
            distances(1, i) = points(i) - 0.5 * (sumLengths(j) + sumLengths(j+1));
        end

    end
        
end

end
