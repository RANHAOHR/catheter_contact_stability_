function [area]=compute_area(direction_angle, area_vec)

    area = sqrt( (area_vec(1) * direction_angle(1))^2 + (area_vec(2) * direction_angle(2))^2 +(area_vec(3) * direction_angle(3))^2 );
end