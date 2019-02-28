function [force, jacobian] = contact_force(obj, jointAngles, currents, externalWrenches)
%
% [force, jacobian] = contact_force(jointAngles, currents, externalWrenches, epsilon)
%
% This function calculates the contact force. It also returns the contact Jacobian
% for your convinience.
%
% Inputs:
% jointAngles is a joint angle vectors.
% currents is an actuator input current vector.
% externalWrenches is an external force.
%
% Outputs:
% force is the contact force
% jacobian is the contact Jacobian (see eqn 5.14 MLS)
%
    tolSurfaceConstraint = 0.5; % obj.tolSurfaceConstraint
    % If the surface is defined and the end-effector is on the surface.
    if (~isempty(obj.surface) && abs(obj.surface.distance(obj.tip_position(jointAngles))) < tolSurfaceConstraint)
        jacobian = contact_jacobian(obj, jointAngles);
        force = pinv(jacobian') * (-obj.stiffnessMatrix * jointAngles + ...
            obj.joint_torques(jointAngles, currents, externalWrenches));
    else
        force = zeros(3, 1);
    end

end

function J = contact_jacobian(obj, q)
%
% J = contact_jacobian(q, param)
%
% Calculate the contact Jacobian (generally know as hand Jacobian).
% See Euation 5.14 in [1] for more details.
%
% Input:
% q is the joint angle vector.
%
% Output:
% J is the Hand Jacobian ().
%

    % Get parameters.
    g_sc = obj.surface.tangent_frame(obj.surface.projection(obj.tip_position(q)));  % Assume tip is on the surface
    g_sf = eye(4);
    p0 = obj.jointPositions;

    % Define wrench basis.
    B_c = [1 0 0; 0 1 0; 0 0 1; zeros(3)];

    % Initialize Jacobian from base to tip frame.
    J_sf = zeros(6, 3 * obj.numJoints);

    % Calculate the finger Jacobian.
    for n = 1:1:obj.numJoints
        w = q(3 * n - 2 : 3 * n);
        p = [0; 0; p0(1, n)];
        g = se3rot(w, p, 1);

        % Calculate the twists analytically.
        [w1, w2, w3] = obj.velocity_axes(w);
        xi1 = twistr(w1, p);
        xi2 = twistr(w2, p);
        xi3 = twistr(w3, p);

        % Calculate the finer Jacobian.
        J_sf(:, 3 * n - 2) = adj(g_sf) * xi1;
        J_sf(:, 3 * n - 1) = adj(g_sf) * xi2;
        J_sf(:, 3 * n) = adj(g_sf)* xi3;
        g_sf = g_sf * g;
    end

    % Assemble the contact (or hand) Jacobian eq5.14 MLS
    J = B_c' * adjinv(g_sc) * J_sf;

end

