function [e_th, e_ph] = get_spherical_basis_vec(dir)
    theta = acos(dir(3)); phi = atan2(dir(2), dir(1));
    [e_th, e_ph] = get_spherical_basis(theta, phi);
end