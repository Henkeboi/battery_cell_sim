% TODO: Derivative of ocp is missing from the expression

phi = 0.5;
ovc_pos = 4.19829 + 0.0565661 * tanh(-14.55460 * phi + 8.60942) - 0.0275479 * ( 1 / ((0.998432 - phi) ^ 0.4924656) - 1.90111) - 0.157123 * exp(-0.04738 * phi ^ 6) + 0.810239 * exp(-40 * (phi - 0.133875)); 
for i = 1 : size(s, 2)
    numerator = L_pos * sqrt(asymmetric_charge_transfer_pos / (resistivity_pos + ionic_conductivity_pos));
    denominator_1 = R_solid_electrolyte_pos;
    denominator_2 = R_pos * ovc_pos / (F * diffusivity_pos);
    denominator_3 = 1.0;
    denominator_4 = R_solid_electrolyte_pos * sqrt(s(i) / diffusivity_pos);
    denominator_5 = (exp(2.0 * R_pos / sqrt(s(i) / diffusivity_pos)) + 1.0) / (exp(2.0 * R_pos / sqrt(s(i) / diffusivity_pos)) - 1.0);
    disp(denominator_4);
    denominator = sqrt(denominator_1 + denominator_2 * (1.0 / (denominator_3 - denominator_4 * denominator_5)));
    nu_pos(i) = numerator / denominator;
    disp(nu_pos(i))
end

% ocv_neg = -0.16 + 1.32 * exp(-3 * phi) + 10 * exp(-2000 * phi);
% ocv_neg_derivative = 10 * 1000 * exp(-2000 * phi);
% disp(ocv_neg_derivative);

