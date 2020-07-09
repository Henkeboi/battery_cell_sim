function [tf_potse, res0, D] = tf_potse(cse, z_coordinates, T_len, sampling_f, electrode, const)
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
        sigma = const.sigma_neg;
        kappa = const.kappa_neg;
        beta = Rs * sqrt(s / Ds);
        eps = const.porosity_solid_neg;
        Uovc_d = calculate_ocv_derivative_neg(cse, const);
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
        sigma = const.sigma_pos;
        kappa = const.kappa_pos;
        beta = Rs * sqrt(s / Ds);
        eps = const.porosity_solid_pos;
        Uovc_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end
    z = z_coordinates(1); % choose z(1)

    % Calculate nu
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uovc_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));

    % Calculate TF potential between solid and electrolyte
    tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z -1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    tf0_potse = (L * z * z * (kappa + sigma) - 2 * L * z * kappa + L * kappa) / (2 * A * kappa * sigma) .* ones(1, size(s, 2)); % Found using maple.
    tf_potse = tf_potse - Uovc_d ./ (eps * A * F * L * s);
    if electrode == 'pos'
        tf_potse = -tf_potse;
    end
    res0 = 0;
    D = 0;

    for i = 1 : size(s, 2)
        if isnan(tf_potse(1, i)) && s(1, i) == 0
            tf_potse(1, i) = tf0_potse(1, i);
        end
    end

    if any(isnan(tf_potse))
        error("NAN in tf_potse");
    end
end
