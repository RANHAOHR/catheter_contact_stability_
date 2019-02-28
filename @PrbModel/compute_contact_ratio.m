function sigma = compute_contact_ratio(obj, f_c_)
    sigma = sqrt(f_c_(1) * f_c_(1) + f_c_(2) * f_c_(2)) / f_c_(3);
end