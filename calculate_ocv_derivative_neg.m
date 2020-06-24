function [ocv_derivative_neg] = calculate_ocv_derivative_neg(cse_neg, const)
    solid_max_c_neg = const.solid_max_c_neg;
    phi = cse_neg / solid_max_c_neg;
    % Derivative of ocv_neg with respect to the surface concentration of lithium in the negative electrolyte with the inital solid concentration substituted.
    ocv_derivative_neg = -1.32 * 3.0 * exp(-3.0 * phi) - 10.0 * 2000.0 * exp(-2000.0 * phi);
end
