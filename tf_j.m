function [tf_j] = phi_se_tf(cse, z_coordinates, const)
    sampling_f = 200;
    T = 1 / sampling_f;
    min_T_len = 4;
    num_samples = 2 ^ (ceil(log2(sampling_f * min_T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    cse_neg = cse;
    z = 0.1; % z coordinate
    Uovc_d = calculate_ocv_derivative_neg(cse_neg, const);

    % Declare params
    Rse = const.R_ct_neg + const.R_film_neg;
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

    % Calculate nu
    nu_num = L * sqrt((alpha / sigma) + (alpha / kappa));
    nu_den = sqrt(Rse + (Rs * Uovc_d) / (F * D)) * sqrt(1 ./ (1 - beta .* coth(beta)));
    nu = nu_num ./ nu_den;

    tf_phi = tf_pot_se(cse_neg, z_coordinates, const);
    tf_j = nu .* tf_phi .* nu / (alpha * F * L * L * ((1 / kappa) + (1 / sigma)));
    tf_j0 = 1 / (L * alpha * F * A); % Found using maple.
    for i = 1 : size(s, 2)
        if isnan(tf_j(1, i)) && s(1, i) == 0
            tf_j(1, i) = tf_j0;
        else if isnan(tf_j(1, i))
            tf_j(1, i) = 0;
        end
    end
end
