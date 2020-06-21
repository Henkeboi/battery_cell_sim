function [nu_neg] = calculate_nu_neg(s, ocv_derivative_neg, sigma_eff_neg, k_eff_neg, L_neg, asymmetric_charge_transfer_neg, resistivity_neg, ionic_conductivity_neg, R_solid_electrolyte_neg, R_neg, diffusivity_neg, F)
    for i = 1 : size(s, 2)
        numerator = L_neg * sqrt(asymmetric_charge_transfer_neg / (resistivity_neg + ionic_conductivity_neg));
        denominator_1 = R_solid_electrolyte_neg;
        denominator_2 = R_neg * ocv_derivative_neg / (F * diffusivity_neg);
        denominator_3 = 1.0;
        denominator_4 = R_solid_electrolyte_neg * sqrt(s(i) / diffusivity_neg);
        denominator_5 = (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) + 1.0) / (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) - 1.0);
        denominator = sqrt(denominator_1 + denominator_2 * (1.0 / (denominator_3 - denominator_4 * denominator_5)));
        nu_neg(i) = numerator / denominator;
    end
end
