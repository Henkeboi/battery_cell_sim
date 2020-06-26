function [ocv_neg] = calculate_ocv_neg(cse_neg, const)
    phi = cse_neg / const.solid_max_c_neg;
    ocv_neg = -0.16 + 1.32 * exp(-3.0 * phi) + 10.0 * exp(-2000.0 * phi);
end
