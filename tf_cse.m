function [tf_cse, res0, D, sampling_f, T_len] = tf_cse(cse, z_coordinates, const, electrode)
    sampling_f = 4000;
    T = 1 / sampling_f;
    T_len = 100;
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
        sigma = const.sigma_pos;
        kappa = const.kappa_pos;
        beta = Rs * sqrt(s / Ds);
        eps = const.porosity_solid_pos;
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    % Calculate nu
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));

    z = 0.0; % Temporarly choose z = 0.0
    % res0 = 1 ./ (eps * A * F * L);
    D = 0;
    % tf_cse = (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) * Rs .* nu ./ (alpha * F * L * A * Ds * (kappa + sigma) * sinh(nu) .* (1 - beta .* coth(beta)));
    % tf_cse0 = (kappa * Rs * (z - 1) + Rs * sigma) / (alpha * F * L * A * (kappa + sigma)); % Found using Maple.


    cse_j = Rs/Ds./(1-beta.*coth(beta));
    res0 = -3/(A*alpha*F*L*Rs);
    tf_cse0 = (-6*Uocv_d*Rs*kappa*sigma +5*alpha*Ds*F*L^2*((2-6*z+3*z.^2)*kappa +(3*z.^2-1)*sigma))./ (30*A*alpha*Ds*Uocv_d*F*L*kappa*sigma);
    tf_cse = nu./(alpha*F*L*A*(kappa+sigma).*sinh(nu)).*(sigma*cosh(nu.*z) + kappa*cosh(nu.*(z-1)));
    tf_cse = tf_cse.*cse_j;    % Convert tf_j to tf_cse
    tf_cse = tf_cse - res0./s; % Remove pole at origin



    % tf_cse = tf_cse + res0 ./ s;
    for i = 1 : size(tf_cse, 2)
        if isnan(tf_cse(1, i)) && s(i, 1) == 0
            tf_cse(1, i) = tf_cse0;
        end
    end

    if electrode == 'pos'
        tf_cse = -tf_cse;
        res0 = -res0;
    end

    if any(isnan(tf_cse))
        error("NAN In tf_cse");
    end
end
