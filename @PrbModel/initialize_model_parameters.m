function initialize_model_parameters(obj, stiffnesses)
%
% initialize_model_parameters(obj, numJoints, linkLengths, stiffnesses)
%
% Initialize parameters exclusively used by the PRB model.
%
% Inputs:
% numJoints is the number of joints.
% linkLengths is the lengths of the links.
% stiffnesses is the stiffnesses of the joints.
%

% Validate stiffness vector if given
if nargin == 2
    assert(size(stiffnesses, 1) == 1);
    assert(size(stiffnesses, 2) == 2 * numJoints);
    
% Calculate stiffness vector if not given
else
    % Moment of inertia of area
    moment_inertia = 0.25 * pi * (obj.outerRadius^4 - obj.innerRadius^4);
    % Polar moment of inertia of area
    polar_moment_inertia = 0.5 * pi * (obj.outerRadius^4 - obj.innerRadius^4);
    % Stiffness vector
    stiffnesses = zeros(1, 2 * obj.numJoints);
    
    for i = 1 : 1 : obj.numJoints
        stiffnesses(1, 2 * i - 1) = obj.youngsModulus * moment_inertia / obj.linkLengths(1, i + 1);
        stiffnesses(1, 2 * i) = obj.shearModulus * polar_moment_inertia / obj.linkLengths(1, i + 1);
    end
    
end

% Store parameters
obj.jointPositions =  obj.linkLengths(1:obj.numJoints) * triu(ones(obj.numJoints));
obj.stiffnessMatrix = diag(kron(eye(obj.numJoints), [1, 0; 1, 0; 0, 1]) * stiffnesses');
obj.linkMasses = obj.massPerLength * obj.linkLengths(1, 2:end);
obj.actuatorLinks = obj.get_links(obj.actuatorPositions);
obj.linkPositions = zeros(1, obj.numJoints);
obj.linkConfigurations = zeros(4, 4, obj.numJoints);

for i = 1 : 1 : obj.numJoints
    obj.linkPositions(1, i) = obj.linkLengths(1) + sum(obj.linkLengths(2:i+1)) - 0.5 * obj.linkLengths(i+1);
    obj.linkConfigurations(1:4, 1:4, i) = [eye(3), [0; 0; obj.linkPositions(1, i)]; 0, 0, 0, 1];
end

end
