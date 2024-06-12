function [sigma_stress] = membrane_bending_stress(epsilon_strain,E_modulus)

    sigma_stress = epsilon_strain.*E_modulus;

end

