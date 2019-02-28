function [g_su, J_su, J_u, B_b, B, P] = actuation_maps(obj, q)
%
% [g_su, J_su, J_u, B_b, B, P] = actuation_maps(q)
%
% Calculates conf of actuation frame, actuator Jacobian, compact form of 
% magnetic field matrix, and input current projection matrix.
%
% Input:
% q is a joint angle vector.
%
% Outputs:
% g_su is the configuration of the actuator with respect to the
%  base frame.
% J_su is the body manipulator Jacobian of the actuator.
% J_u is the bottom half of the manipulator Jacobian.
% B_b is MRI's magnetic field in the body frame.
% B maps input current in the plane (R^2) orthogonal to the
%  magnetic field to actuator torque vector in body frame.
% P maps the input current from R^3 to the input in R^2.
%

% obj.parameterseters
p0 = obj.jointPositions;
w0 = obj.jointDirections;
B_s = obj.magneticField;  % Magnetic field in the spatial frame
NA = obj.turnAreaMatrices;
CC = obj.coilAlignmentMatrices;
n = obj.numJoints;
m = obj.numActuators;
% Initialize matrices
actuatorTorqueDofs = 2;
g_su = obj.actuatorConfigurations;  % Initial conf
J_su = zeros(6, obj.jointDofs * n, m);  % Init body Jacobian
J_u = zeros(obj.actuatorDofs * m, obj.jointDofs * n);  % Init rotational component of body Jacobian
B = zeros(obj.actuatorDofs * m, actuatorTorqueDofs * m);  % Init compact representation of input matrix
P = zeros(actuatorTorqueDofs * m, obj.actuatorDofs * m);  % Init mapping from compact to non-compact input
actuator_links = obj.actuatorLinks;
% Calculate body Jacobians
for i = 1:1:obj.numActuators
    % The joint immediately below the link
    joint_number = actuator_links(1, i);
    for j = joint_number:-1:1
        w = q(3*j-2:3*j);
        p = [0; 0; p0(1, j)];
        g_su(:, :, i) = se3rot(w, p, 1) * g_su(:, :, i);
        [w_x, w_y, w_z] = obj.velocity_axes(w);
        xi_x = twistr(w_x, p);
        xi_y = twistr(w_y, p);
        xi_z = twistr(w_z, p);
        J_su(:, 3*j-2, i) = adjinv(g_su(:, :, i)) * xi_x;
        J_su(:, 3*j-1, i) = adjinv(g_su(:, :, i)) * xi_y;
        J_su(:, 3*j, i) = adjinv(g_su(:, :, i)) * xi_z;
    end
    % Calculate compact form of magnetic field cross product
    B_b = g_su(1:3, 1:3, i)' * B_s;           % Magnetic field in the body frame
    % Take out null space of magnetic field via SVD
    try
        [U, S, V] = svd(-vector2skewsym(B_b) * CC(:, :, i) * NA(:, :, i)); % SVD of -B_b * NA (minus sign from reversing order of cross product)
    catch
        warning('Problem with SVD of input matrix');
    end
    % Input matrices
    B(3*i-2:3*i, 2*i-1:2*i) = U(1:3, 1:2) * S(1:2, 1:2);            % Calculate the compact representation of B_b
    P(2*i-1:2*i, 3*i-2:3*i) = V(1:3, 1:2)';                      % Input projection
    % Actuator Jacobian
    J_u(3*i-2:3*i, :) = J_su(4:6, :, i);                     % Bottom part of the Jacobian
end

end
