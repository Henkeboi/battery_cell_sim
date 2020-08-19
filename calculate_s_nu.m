function [s, nu] = calculate_nu(cse, T_len, sampling_f, electrode, const)
    T = 1 / sampling_f;
    num_samples = 2 ^ (ceil(log2(sampling_f * T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    if electrode == 'neg'
        Rct = const.R_ct_neg;
        Rfilm = const.R_film_neg;
        Rse = Rct + Rfilm;
        Rs = const.radius_neg;
        L = const.L_neg;
        Ds = const.diffusivity_neg;
        alpha = const.alpha_neg;
        kappa = const.kappa_eff_neg;
        sigma = const.sigma_eff_neg;
        as = const.as_neg;
        A = const.A_neg;
        F = const.F;
        beta = Rs * sqrt(s / Ds);
        Uocv_d = calculate_ocv_derivative_neg(cse, const);
    elseif electrode == 'pos'
        Rct = const.R_ct_pos;
        Rfilm = const.R_film_pos;
        Rse = Rct + Rfilm;
        Rs = const.radius_pos;
        L = const.L_pos;
        Ds = const.diffusivity_pos;
        alpha = const.alpha_pos;
        kappa = const.kappa_eff_pos;
        sigma = const.sigma_eff_pos;
        as = const.as_pos;
        A = const.A_pos;
        F = const.F;
        beta = Rs * sqrt(s / Ds);
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end
    nu = L * sqrt((as / sigma + as / kappa) ./ (Rse + Rs * Uocv_d ./ (F * Ds - F * Ds * beta .* coth(beta))));
    nu(1) = 0.0;

    if any(isnan(nu))
        error("NAN in calculate nu")
    end
end
