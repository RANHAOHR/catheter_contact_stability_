function jointTorques = joint_torques(obj, jointAngles, currents, externalWrenches)
%
% jointTorques = joint_torques(jointAngles, currents, externalWrenches)
%
% This function calculates the joint torque vector due to actuation 
% and external wrenches. The joint torque vector is calculated by 
% propagating torques down the body of the catheter. While this method 
% is messier, it is faster than calculating the body Jacobian first 
% then multiplying it with the wrenches.
%
% Inputs:
% jointAngles is a joint angle vector.
% currents is a current vector.
% externalWrenches is an array of external wrenches acting on the
%  center of mass of the links. Its dimension is 6-by-N, where 6 is
%  the dimension of each wrench, and N is the number of joints.
%
% Output:
% jointTorques is the resulting joint torque vector.
%

% Initialize joint torques
jointTorques = zeros(obj.get_dofs(), 1);

% Temporary joint torques used to propagate joint torques down the links.
jointWrench = zeros(6, 1);

% Calculate configurations of points of interest.
jointConfigurations = obj.point_configurations(jointAngles, obj.jointPositions);
linkConfigurations = obj.point_configurations(jointAngles, obj.linkPositions);
actuatorConfigurations = obj.point_configurations(jointAngles, obj.actuatorPositions);

% Actuator index.
actuatorIndex = obj.numActuators;

for jointIndex = obj.numJoints : -1 : 1
    
    % The ith link is above the ith joint.
    linkIndex = jointIndex;
    
    % Indices associated with this joint.
    jointIndices = obj.jointDofs * (jointIndex - 1) + 1 : obj.jointDofs * jointIndex;
    
    % If this is not the last joint, add torque from previous joint.
    if jointIndex < obj.numJoints
        % Configuration and adjoint matrices of center of mass with respect to the joint.
        configurationDA = jointConfigurations(:, :, jointIndex + 1) \ jointConfigurations(:, :, jointIndex);
        adjointDA = adj(configurationDA);
        jointWrench = adjointDA' * jointWrench;
    end
    
    % If this link contains an actuator, calculate it joint torques.
    if actuatorIndex > 0 && linkIndex == obj.actuatorLinks(actuatorIndex)
        
        % Currents going though this actuator
        actuatorCurrents = currents(obj.actuatorDofs * (actuatorIndex - 1) + 1 : obj.actuatorDofs * actuatorIndex);
        
        % Configuration and adjoint matrices of the actuator.
        configurationCA = actuatorConfigurations(:, :, actuatorIndex) \ jointConfigurations(:, :, jointIndex);
        adjointCA = adj(configurationCA);
        
        % Update joint torques.
        magneticMomentsBody = obj.coilAlignmentMatrices(:, :, actuatorIndex) * ...
            obj.turnAreaMatrices(:, :, actuatorIndex) * actuatorCurrents;
        magneticFieldBody = actuatorConfigurations(1:3, 1:3, actuatorIndex)' * ...
            obj.magneticField;
        actuatorWeightBody = actuatorConfigurations(1:3, 1:3, actuatorIndex)' * ...
            obj.gravity * obj.coilMasses(1, actuatorIndex);
        jointWrench = jointWrench + ...
            adjointCA' * [actuatorWeightBody; cross(magneticMomentsBody, magneticFieldBody)];
        
        % Update actuatorIndex.
        actuatorIndex = actuatorIndex - 1;
    end
    
    % Configuration and adjoint matrices of center of mass with respect to the joint.
    configurationBA =  linkConfigurations(:, :, linkIndex) \ jointConfigurations(:, :, jointIndex);
    adjointBA = adj(configurationBA);
    
    % Calculate joint torques due to gravity and external wrenches.
    linkWeightWrenchBody = [linkConfigurations(1:3, 1:3, linkIndex)' * obj.gravity * obj.linkMasses(1, linkIndex); zeros(3, 1)];
    externalWrenchBody = adj([linkConfigurations(1:3, 1:3, linkIndex), zeros(3, 1); 0, 0, 0, 1])' * ...
        externalWrenches(:, jointIndex);
    jointWrench = jointWrench + adjointBA' * (linkWeightWrenchBody + externalWrenchBody);
    
    % Calculate joint torques.
    [twistX, twistY, twistZ] = obj.velocity_axes(jointAngles(jointIndices));
    jointTorques(jointIndices, 1) = [twistX, twistY, twistZ]' * jointWrench(4:6, 1);
    
end

end
