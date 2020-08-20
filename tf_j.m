function [tf_j, res0, D] = tf_j(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);

    % Declare params
    if electrode == 'neg'
        Rct = const.R_ct_neg;
        Rfilm = const.R_film_neg;
        Rse = Rct + Rfilm;
        Rs = const.radius_neg;
        as = const.as_neg;
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
        as = const.as_pos;
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

    z = 0.4;
    tf_j = nu .* (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (as * A * F * L * (kappa + sigma) * sinh(nu));
    tf_j(1) = 1 / (as * F * L * A);
    [row, col] = find(isnan(tf_j));
    tf_j(row, col) = 0;
    res0 = 0;
    D = NaN;

    if electrode == 'pos'
        tf_j = -tf_j;
    end
   
    if any(isnan(tf_j))
        error("NAN in tf_j");
    end
end
