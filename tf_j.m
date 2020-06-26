function [tf_j, res0, D, sampling_f, T_len] = phi_se_tf(cse, z_coordinates, const, electrode)
    sampling_f = 400;
    T = 1 / sampling_f;
    T_len = 10;
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
        Rct = const.R_ct_pos;
        Rfilm = const.R_film_pos;
        Rse = Rct + Rfilm;
        Rs = const.radius_pos;
        L = const.L_pos;
        F = const.F;
        D = const.diffusivity_pos;
        A = const.A_pos;
        alpha = const.alpha_pos;
        sigma = const.sigma_pos;
        kappa = const.kappa_pos;
        beta = Rs * sqrt(s / D);
        eps = const.porosity_solid_pos;
        Uovc_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    % Calculate nu
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uovc_d / F / D * (1 ./ (Rs * beta .* coth(beta))));

    % Calculate TF potential between solid and electrolyte.
    tf_phi = zeros(size(s, 2), size(z_coordinates, 2));
    for i = 1 : size(z_coordinates, 2)
        z = z_coordinates(1, i); 
        if electrode == 'neg'
            tf_phi = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu)); % PHI / Iapp.
            tf_j(:, i) = nu .* tf_phi .* nu / (alpha * F * L * L * (1 / kappa + 1 / sigma));
        else
            tf_phi = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu)); % PHI / Iapp.
            tf_j(:, i) = -nu .* tf_phi .* nu / (alpha * F * L * L * (1 / kappa + 1 / sigma));
        end
    end

    if electrode == 'neg'
        tf_j0 = 1 / (L * alpha * F * A); % Found using Maple.
    else
        tf_j0 = -1 / (L * alpha * F * A); % Found using Maple.
    end
    for i = 1 : size(z_coordinates, 2) 
        for j = 1 : size(s, 2)
            if isnan(tf_j(j, i)) && s(1, j) == 0
                tf_j(j, i) = tf_j0;
            else if isnan(tf_j(j, i))
                tf_j(j, i) = 0;
            end
        end
    end

    % Find res0
    res0 = zeros(size(z_coordinates, 2), size(z_coordinates, 2));

    % Find the unit impulse response where t = 0 and s -> inf.
    nu_inf = L * sqrt((alpha / kappa + alpha / sigma) / (Rct + Rfilm));  % From Gregory Plett.
    D = zeros(size(nu_inf, 2));
    
    if any(isnan(tf_j))
        error("NAN in tf_j");
    end
end
