function [F_e_b, bodyJacobian] = external_drag_(obj, jointAngles, externalWrenches)
%
% [jointTorques, bodyJacobian] = external_joint_torque(jointAngles, externalWrenches)
%
% Calculate joint torques due to gravity.
%
% Inputs:
% jointAngles is joint angle (column) vector.
% externalWrenches is either disturbance wrench in R^6 for each disturbance
%  point.
%
% Outputs:
% jointTorques is joint torque (column) vector.
% bodyJacobian is a three-dimensional array where the body manipulator 
%  Jacobians are stacked along the third dimension. It is indexed by 
%  J(row, column, jointNumber).
%

% Initialize containers.
F_e_b = zeros(obj.jointDofs * obj.numJoints, 1);
bodyJacobian = zeros(6, obj.jointDofs * obj.numJoints, obj.numJoints);

% Calculate body Jacobians and joint torques,
for jointNumber = 1:1:obj.numJoints
    % The corresponding link.
    linkNumber = jointNumber;
    % Initial configuration of center of mass of the link.
    linkConfiguration = obj.linkConfigurations(:, :, linkNumber);
    
    % Calculate body Jacobian.
    for i = jointNumber:-1:1
        iJointAngles = jointAngles(3 * i - 2 : 3 * i, 1);
        jointPosition = [0; 0; obj.jointPositions(1, i)];
        linkConfiguration = se3rot(iJointAngles, jointPosition, 1) * linkConfiguration;
        [axisX, axisY, axisZ] = obj.velocity_axes(iJointAngles);
        twistX = twistr(axisX, jointPosition);
        twistY = twistr(axisY, jointPosition);
        twistZ = twistr(axisZ, jointPosition);
        bodyJacobian(:, 3 * i - 2, jointNumber) =  adjinv(linkConfiguration) * twistX;
        bodyJacobian(:, 3 * i - 1, jointNumber) =  adjinv(linkConfiguration) * twistY;
        bodyJacobian(:, 3 * i, jointNumber) =  adjinv(linkConfiguration) * twistZ;
    end
    
%     % External wrench in body frame.
%     externalWrenchBody = [linkConfiguration(1:3, 1:3)', zeros(3, 3); zeros(3, 3), linkConfiguration(1:3, 1:3)'] * externalWrenches(:, jointNumber);
% 
%     % Add to joint torques.
%     F_e_b = F_e_b + bodyJacobian(:, :, jointNumber)' * externalWrenchBody;

    externalWrenchBody = [linkConfiguration(1:3, 1:3)', zeros(3, 3); zeros(3, 3), linkConfiguration(1:3, 1:3)'] * externalWrenches(:,:,jointNumber);

    % Add to joint torques.
    F_e_b = F_e_b + bodyJacobian(:, :, jointNumber)' * externalWrenchBody;

end

end
