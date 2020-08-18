function [tf_cse, res0, D] = tf_cse(cse, z_coordinates, T_len, sampling_f, electrode, const)
    T = 1 / sampling_f;
    num_samples = 2 ^ (ceil(log2(sampling_f * T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    % Declare params
    if electrode == 'neg'
        Rct = const.R_ct_neg;
        Rfilm = const.R_film_neg;
        Rse = Rct + Rfilm;
        Rs = const.radius_neg;
        L = const.L_neg;
        F = const.F;
        Ds = const.diffusivity_neg;
        A = const.A_neg;
        alpha = const.alpha_neg;
        sigma = const.sigma_eff_neg;
        kappa = const.kappa_eff_neg;
        beta = Rs * sqrt(s / Ds);
        eps = const.eps_s_neg;
        Uocv_d = calculate_ocv_derivative_neg(cse, const);
    elseif electrode == 'pos'
        Rct = const.R_ct_pos;
        Rfilm = const.R_film_pos;
        Rse = Rct + Rfilm;
        Rs = const.radius_pos;
        L = const.L_pos;
        F = const.F;
        Ds = const.diffusivity_pos;
        A = const.A_pos;
        alpha = const.alpha_pos;
        sigma = const.sigma_eff_pos;
        kappa = const.kappa_eff_pos;
        beta = Rs * sqrt(s / Ds);
        eps = const.eps_s_pos;
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));
    z = 0.5;
    tf_cse = (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) * Rs .* nu;
    tf_cse0 = 0;
    tf_cse(1) = tf_cse0;

    res0 = 0;
    D = 0;

    if electrode == 'pos'
        tf_cse = -tf_cse;
    end

    if any(isnan(tf_cse))
        error("NAN In tf_cse");
    end
end
