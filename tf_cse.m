function [tf_cse, res0, D, sampling_f, T_len] = phi_se_tf(cse, z_coordinates, const, electrode)
    sampling_f = 200;
    T = 1 / sampling_f;
    T_len = 1;
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
        D = const.diffusivity_neg;
        A = const.A_neg;
        alpha = const.alpha_neg;
        sigma = const.sigma_neg;
        kappa = const.kappa_neg;
        beta = Rs * sqrt(s / D);
        eps = const.porosity_solid_neg;
        Uovc_d = calculate_ocv_derivative_neg(cse, const);
    elseif electrode == 'pos'
        Rct = const.R_ct_neg;
        Rfilm = const.R_film_neg;
        Rse = Rct + Rfilm;
        Rs = const.radius_neg;
        L = const.L_neg;
        F = const.F;
        D = const.diffusivity_neg;
        A = const.A_neg;
        alpha = const.alpha_neg;
        sigma = const.sigma_neg;
        kappa = const.kappa_neg;
        beta = Rs * sqrt(s / D);
        eps = const.porosity_solid_neg;
        Uovc_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end
    
    % Calculate nu
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uovc_d / F / D * (1 ./ (Rs * beta .* coth(beta))));
    
    % Calcualte lithium concentration on the electrode surface TF.
    z = 0.5; % Temporary choosen z-location
    tf_cse = (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) .* (Rs * nu) ./ (alpha * F * L * A * D * (kappa + sigma) * sinh(nu) .* (1 - sqrt(beta) .* coth(beta)));
    tf_cse0 = (kappa * Rs * (z - 1) + Rs * sigma) / (alpha * F * L * A * (kappa + sigma)); % Found using Maple.
    res0 = 0;
    D = 0;
    
    for i = 1 : size(s, 2)
        if isnan(tf_cse(1, i)) && s(1, i) == 0
            tf_cse(1, i) = tf_cse0;
        end
    end

    if any(isnan(tf_cse))
        error("NAN In tf_cse")
    end
end
