% TODO: Derivative of ocp is missing from the expression
for i = 1 : size(s, 2)
    numerator =  L_neg * sqrt(asymmetric_charge_transfer_neg / (resistivity_neg + ionic_conductivity_neg));
    denominator = sqrt(R_solid_electrolyte_neg + R_neg / (F * diffusivity_neg) * (1 / (1 - R_solid_electrolyte_neg * coth(R_neg / sqrt(s(i) / diffusivity_neg)) * sqrt(s(i) / diffusivity_neg))));
    nu_neg(i) = numerator / denominator;
end
