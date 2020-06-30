% v(t) = potpos(0, t).
% v is needed to bias the positive electrode potential
% x-coordiante = 0 at neg electrode current collector and L at pos electrode current collector.
% z-coordinate = 0 at the current collector and z = 1 at the electrolyte interface.
% U in z coordinates eg z = x / L.
% voltage(t) = Ø_pos(0, t) - Ø_neg(0, t) in z coordinates.
% We need to know V to debias debiased_Ø_neg.
% Linearization of the Butler-Volmer equiation: 

% TODO: Currently using the initial or wrong values for concentration, j and phi.

function [v] = calculate_voltage(params, const)
    R = const.R;
    F = const.F;
    temp = const.temp;
    j0_pos = const.j0_pos;
    j0_neg = const.j0_neg;
    k_neg = const.reaction_rate_neg;
    k_pos = const.reaction_rate_pos;
    ce_neg = const.initial_electrolyte_concentration_neg;
    ce_pos = const.initial_electrolyte_concentration_pos;
    cs_max_neg = const.solid_max_c_neg;
    cs_max_pos = const.solid_max_c_pos;
    Rf_neg = const.R_film_neg;
    Rf_pos = const.R_film_pos;
    Rse_neg = const.R_solid_electrolyte_neg;
    Rse_pos = const.R_solid_electrolyte_pos;
    Uocp_neg = calculate_ocv_neg(params.cse_neg, const); % Calculate cse_neg at z = 0
    Uocp_pos = calculate_ocv_pos(params.cse_pos, const); % Calculate cse_pos at z = 0
    beta_pos = sqrt(ce_pos * (cs_max_pos - params.cse_pos) * params.cse_pos);
    beta_neg = sqrt(ce_neg * (cs_max_neg - params.cse_neg) * params.cse_neg);

    n_pos = 2 * R * temp / F * asinh(params.j_pos / (2 * k_pos * beta_pos)); % Calculate at z = 0 
    n_neg = 2 * R * temp / F * asinh(params.j_neg / (2 * k_neg * beta_neg)); % Calculate at z = 0

    v = F * (Rf_pos * params.j_pos - Rf_neg * params.j_neg) + params.pote1 + params.pote2 + Uocp_pos - Uocp_neg
end
