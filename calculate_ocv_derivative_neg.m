cse_neg = 4000.0;  % TODO: Calculate cse_neg.
inital_solid_c_neg = 4000.0; % TODO: Use SOC to estimate.
phi = cse_neg / solid_max_concentration_neg;

% Derivative of ocv_neg with respect to the surface concentration of lithium in the negative electrolyte with the inital solid concentration substituted.
ocv_derivative_neg = -.16 / solid_max_concentration_neg  - 1.32 * 3.0 * exp(-3.0 * phi) / solid_max_concentration_neg - 10.0 * 2000.0 * exp(-2000.0 * phi) / solid_max_concentration_neg
