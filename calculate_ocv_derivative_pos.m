function [ocv_derivative_pos] = calculate_ocv_derivative_neg(cse_pos, const)
    solid_max_c_pos = const.solid_max_c_pos;
    phi = cse_pos / solid_max_c_pos;
    % Derivative of ocv_pos with respect to the surface concentration of lithium in the positive electrolyte with the inital solid concentration substituted.
    ocv_derivative_pos = -0.0565661 * 14.5546 * (sech(-14.5546 * phi + 8.60942)) ^ 2 - 0.01356 / ((0.998432 - phi) ^ 1.4924656) + 0.157123 * exp(-0.04738 * phi ^ 6) * 0.04738 * 6 * (phi ^ 5) - 40 * 0.810239 * exp(-40 * (phi - 0.133875))
end
