function [v] = calculate_voltage(j_neg, j_pos, cse_neg, cse_pos, ce_neg, ce_pos, pote_1, pote_2, const)
    F = const.F;
    R = const.R;
    T = const.temp;
    kappa_norm_neg = const.kappa_norm_neg; % k0_neg
    kappa_norm_pos = const.kappa_norm_pos; % k0_pos
    ce0_neg = const.ce0_neg;
    ce0_pos = const.ce0_pos;
    cs_max_neg = const.solid_max_c_neg;
    cs_max_pos = const.solid_max_c_pos;
    R_film_neg = const.R_film_neg;
    R_film_pos = const.R_film_pos;

    % This assumes alpha equal to 0.5
    n_neg = 2 * R * T / F * asinh(j_neg / (2 * kappa_norm_neg * sqrt(ce_neg * (cs_max_neg - cse_neg) * cse_neg)));
    n_pos = 2 * R * T / F * asinh(j_pos / (2 * kappa_norm_pos * sqrt(ce_pos * (cs_max_pos - cse_pos) * cse_pos)));
    
    v = F * (R_film_pos * j_pos - R_film_neg * j_neg) + pote_1 + n_pos - n_neg + pote_2 + calculate_ocv_pos(cse_pos, const) - calculate_ocv_neg(cse_neg, const);
end
