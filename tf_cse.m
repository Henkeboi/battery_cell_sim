function [tf_cse, res0, D] = tf_cse(cse, z_coordinates, T_len, sampling_f, electrode, const)
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
        cs0 = const.cs0_neg;
        Uocv_d = calculate_ocv_derivative_neg(cs0, const);
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
        cs0 = const.cs0_pos;
        Uocv_d = calculate_ocv_derivative_pos(cs0, const);
    else
        error("Bad electrode selection");
    end

    z = 0.5;
    nu = L * sqrt(as / sigma + as / kappa) ./ sqrt((Rse + Rs * Uocv_d / (F * Ds)) ./ (1 - beta .* coth(beta)));
    %[M, index] = max(nu);
    %nu_inf = L * sqrt((as / sigma + as / kappa) / Rse);
    %nu(index) = nu_inf;

    %tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    %tf_potse0 = L * (z * z * (kappa + sigma) / 2 + kappa * (1 / 2 - z));
    %tf_potse(1) = tf_potse0;
    %tf_j = nu .* nu .* tf_potse / (as * F * L * L * (1 / kappa + 1 / sigma));
    %tf_j(1) = 0.0;
    %tf_cse = Rs / Ds ./ (1 - beta .* coth(beta)) .* tf_j;
    %tf_cse(1) = Rs / Ds * tf_j(1);
    tf_cse = Rs * sigma * nu .* coth(nu * z) ./ (as * F * L * A * Ds * (kappa + sigma) * (1 - beta .* coth(beta)));
    tf_cse = tf_cse + Rs * kappa * nu .* (exp(z) ./ (exp(nu) - 1) + exp(-z * nu) ./ (1 - exp(-nu))) ./ (as * F * L * A * Ds * (kappa + sigma) * (1 - beta .* coth(beta)));
    
    tf_cse(1) = 0;

     

    res0 = 0;
    D = 0;

    if electrode == 'pos'
        tf_cse = -tf_cse;
        res0 = -res0;
    end
    if any(isnan(tf_cse))
        for i = 1 : size(tf_cse, 2)
            if isnan(tf_cse(1, i))
                disp(i)
            end
        end
        error("NAN In tf_cse");
    end
end
