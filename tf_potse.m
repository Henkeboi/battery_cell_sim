function [tf_potse, res0, D] = tf_j(cse, z_coordinates, T_len, sampling_f, electrode, const)
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
        eps = const.porosity_solid_neg;
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
        eps = const.porosity_solid_pos;
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    % Calculate nu
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));
    % Calculate TF potential between solid and electrolyte.
    z = 0;
    tf_potse0 = -Rs / 5 / Ds; % 1 / (L * alpha * F * A); % Found using Maple.
    res0 = (Uocv_d * 3 / Rs) ./ s; % From DRA example 3.
    res0(1, 1) = -Uocv_d * Rs / (5 * Ds); % Analytic solution DRA example 3.
    tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu)) + res0;
    
    for i = 1 : size(tf_potse, 2)
        if isnan(tf_potse(1, i))
            tf_potse(1, i) = tf_potse0;
        end
    end

    res0 = 0;
    D = 0;
    if electrode == 'pos'
        tf_potse = -tf_potse;
        res0 = -res0;
    end
   
    if any(isnan(tf_potse))
        error("NAN in tf_potse");
    end
end
