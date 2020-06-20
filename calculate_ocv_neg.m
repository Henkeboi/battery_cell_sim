cse_neg = 4000; % TODO: calculate cse_neg
phi = cse_neg / solid_max_concentration_neg;
ocv_neg = -0.16 * phi + 1.32 * exp(-3.0 * phi) + 10.0 * exp(-2000 * phi);
