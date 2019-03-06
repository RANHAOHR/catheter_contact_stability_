function [wx, wy, wz] = angle_axes(w_vec)  % From Tipakron's version
% [wx, wy, wz] = velocity_axes(q(1), q(2), q(3))
%
% Calculate velocity axes, wx and wy, of the catheter given its joint
% angles.
%
% Input:
% q is a joint angle vector
%
% Outputs:
% wx is velocity axis in the x-direction
% wy is velocity axis in the y-direction
% wz is velocity axis in the z-direction
%

% Check size
assert(size(w_vec, 1) == 3);
assert(size(w_vec, 2) == 1);

% Store variables for readability
phi = norm(w_vec);
x = w_vec(1, 1);
y = w_vec(2, 1);
z = w_vec(3, 1);

% Initialize veclocity axes
wx = [1; 0; 0];
wy = [0; 1; 0];
wz = [0; 0; 1];

if phi > 1e-16
    wx(1, 1) = (x^4 + x^2 * y^2 + x^2 * z^2 + y^2 * phi * sin(phi) + z^2 * phi * sin(phi)) / phi^4;
    wx(2, 1) = x * y / phi^2 - x * y * sin(phi) / phi^3 - z * cos(phi) / phi^2 + z / phi^2;
    wx(3, 1) = (x^3 * z + x^2 * y * cos(phi) - x^2 * y + x * y^2 * z + x * z^3 - x * z * phi * sin(phi) + y^3 * cos(phi) - y^3 + y * z^2 * cos(phi) - y * z^2) / phi^4;
    wy(1, 1) = (x^3 * y + x^2 * z * cos(phi) - x^2 * z + x * y^3 + x * y * z^2 - x * y * phi * sin(phi) + y^2 * z * cos(phi) - y^2 * z + z^3 * cos(phi) - z^3) / phi^4;
    wy(2, 1) = (x^2 * y^2 + x^2 * phi * sin(phi) + y^4 + y^2 * z^2 + z^2 * phi * sin(phi)) / phi^4;
    wy(3, 1) = -x * cos(phi) / phi^2 + x / phi^2 + y * z / phi^2 - y * z * sin(phi) / phi^3;
    wz(1, 1) = x * z / phi^2 - x * z * sin(phi) / phi^3 - y * cos(phi) / phi^2 + y / phi^2;
    wz(2, 1) = (x^3 * cos(phi) - x^3 + x^2 * y * z + x * y^2 * cos(phi) - x * y^2 + x * z^2 * cos(phi) - x * z^2 + y^3 * z + y * z^3 - y * z * phi * sin(phi)) / phi^4;
    wz(3, 1) = (x^2 * z^2 + x^2 * phi * sin(phi) + y^2 * z^2 + y^2 * phi * sin(phi) + z^4) / phi^4;
end

end
