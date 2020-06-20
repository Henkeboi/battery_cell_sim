cse_pos = 8000;  % TODO: Calcultate solid_concentration_pos
phi = cse_pos / solid_max_concentration_pos;
ocv_pos = 4.192829 - 0.0565661 * tanh(-14.5546 * phi + 8.60942) - 0.0275479 * (1 / ((0.998432 -phi) ^ 0.4924656) - 1.90111) - 0.157123 * exp(-0.04738 * (phi ^ 6)) + 0.810239 * exp(-40 * (phi - 0.133875));
