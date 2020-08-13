function [tf_j, res0] = tf_j(cse, z_coordinates, T_len, sampling_f, electrode, const)
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

    % % Calculate nu
    % nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));

    % % Calculate TF potential between solid and electrolyte.
    % tf_phi = zeros(size(s, 2), size(z_coordinates, 2));
    % tf_j0 = 1 / (L * alpha * F * A); % Found using Maple.
    % for i = 1 : size(z_coordinates, 2)
    %     z = z_coordinates(1, i); 
    %     tf_phi = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu)); % PHI / Iapp.
    %     tf_j(:, i) = nu .* tf_phi .* nu / (alpha * F * L * L * (1 / kappa + 1 / sigma));
    % end

    % for i = 1 : size(z_coordinates, 2) 
    %     for j = 1 : size(s, 2)
    %         if isnan(tf_j(j, i)) && s(1, j) == 0
    %             tf_j(j, i) = tf_j0;
    %         end
    %     end
    % end
    % res0 = 0;


    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));
    z = 0.5;
    tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    tf_potse0 = L * (z * z * (kappa + sigma) / 2 + kappa * (1 / 2 - z));
    tf_potse(1, 1) = tf_potse0;
    tf_j = nu .* nu .* tf_potse / (alpha * F * L * L * (1 / kappa + 1 / sigma));
    tf_j(1, 1) = 0.0;
    res0 = 0;

    if electrode == 'pos'
        tf_j = -tf_j;
        res0 = -res0;
    end
   
    if any(isnan(tf_j))
        error("NAN in tf_j");
    end
end
