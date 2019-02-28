function f = potential_energy(obj, jointAngles, currents, externalWrenches, initialJointAngles)
    
    % Initialize objective function value with internal stiffness potential
    f = 0.5 * (jointAngles - initialJointAngles)' * obj.stiffnessMatrix * (jointAngles - initialJointAngles);
    
    % Add potential energy from coils
    actuatorConfigurations = point_configurations(obj, jointAngles, obj.actuatorPositions);
    
    for i = 1 : 1 : obj.numActuators
        % Magnetic moment aligned with spatial frame
        magneticMoment = actuatorConfigurations(1:3, 1:3, i) * obj.coilAlignmentMatrices(:, :, i) * obj.turnAreaMatrices(:, :, i) * ...
            currents(3 * (i - 1) + 1 : 3 * (i - 1) + 3);
        % Update potential enerfy
        f = f - magneticMoment' * obj.magneticField - ...
            obj.coilMasses(i) * obj.gravity' * actuatorConfigurations(1:3, 4, i);
    end
    
    % Add potential energy from link weights
    linkConfigurations = point_configurations(obj, jointAngles, obj.linkPositions);
    
    for i = 1 : 1 : obj.numJoints
        % Update potential energy
        f = f - obj.linkMasses(1, i) * obj.gravity' * linkConfigurations(1:3, 4, i);
    end
end