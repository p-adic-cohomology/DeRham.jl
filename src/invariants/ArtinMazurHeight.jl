"""
    artin_mazur_height
    computes the Artin-Mazur height of a quartic K3 surface given the coefficients of its L-polynomial
"""
function artin_mazur_height(coeffs, q)
    @assert length(coeffs) == 22
    raw = _normalize_raw_coefficients(coeffs, q)
    boatshape_coeffs = normalized_tate_twist(raw, q, 1)
    for i = 2:11
        if mod(boatshape_coeffs[i], q) != 0
            return i - 1
        end
    end
    return Inf
end
