function [T, g_su, J_su, J_u, B_b, B, P] = actuator_joint_torque(obj, q, u)
%
% [T, g_su, J_su, J_u, B_b, B, P] = actuator_joint_torque(q, u)
%
% Calculate joint torques for a given configuration and control.
%
% Inputs:
% q is a joint angle vector.
% u is an input current.
%
% Outputs:
% T is the joint torques.
% g_su is the configuration of the actuator with respect to the
% base frame.
% J_su is the body manipulator Jacobian of the actuator.
% J_u is the bottom half of the manipulator Jacobian.
% B_b is MRI's magnetic field in the body frame.
% B maps input current in the plane (R^2) orthogonal to the
% magnetic field to actuator torque vector.
% P maps the input current from R^3 to the input in R^2.
%

% Initial gravity vector in spatial frame
g = obj.gravity;
% Initialize joint torques from actuation and coil masses
T = zeros(3 * obj.numJoints, 1);
% Calculate actuation map
[g_su, J_su, J_u, B_b, B, P] = obj.actuation_maps(q);    

% Actuation joint torques
if size(u, 1) == 3 * obj.numActuators
    T = J_u' * B * P * u;
elseif size(u, 1) == 2 * obj.numActuators
    T = J_u' * B * u;    
end

% Gravity joint torques
for m = 1:1:obj.numActuators
    % Gravity in the body frame
    g_b = g_su(1:3, 1:3, m)' * g;
    % Joint torques
    T = T + J_su(:, :, m)'*[obj.coilMasses(1, m) * g_b; 0; 0; 0];    
end

end
