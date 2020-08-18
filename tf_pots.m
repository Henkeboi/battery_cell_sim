function [tf_pots, res0, D] = tf_pots(cse, z_coordinates, T_len, sampling_f, electrode, const)
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
        Ds = const.diffusivity_neg;
        alpha = const.alpha_neg;
        kappa = const.kappa_eff_neg;
        sigma = const.sigma_eff_neg;
        A = const.A_neg;
        F = const.F;
        beta = Rs * sqrt(s / Ds);
        Uocv_d = calculate_ocv_derivative_neg(cse, const);
    elseif electrode == 'pos'
        Rct = const.R_ct_pos;
        Rfilm = const.R_film_pos;
        Rse = Rct + Rfilm;
        Rs = const.radius_pos;
        L = const.L_pos;
        Ds = const.diffusivity_pos;
        alpha = const.alpha_pos;
        kappa = const.kappa_eff_pos;
        sigma = const.sigma_eff_pos;
        A = const.A_pos;
        F = const.F;
        beta = Rs * sqrt(s / Ds);
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    %z = 0.0000000000000001;
    z = 0.5;
    nu = calculate_nu(cse, T_len, sampling_f, electrode, const);
    %nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uocv_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));

    %tf_pots = -L * ((kappa * (cosh(nu) - cosh(nu .* (z - 1)))) + sigma * (1 - cosh(z .* nu) + z .* nu .* sinh(nu))) ./ (A * sigma * (kappa + sigma) .* nu .* sinh(nu));
    %tf_pots0 = -L * z * (1 - z / 2) / (A * sigma);

    tf_pots = -L*(kappa*(cosh(nu) - cosh((z-1).*nu)) +...
        sigma*(1-cosh(z.*nu) + z.*nu.*sinh(nu))) ...
        ./(A*sigma*(kappa + sigma)*...
        nu.*sinh(nu));

    tf_pots0 = L*(z-2).*z/(2*A*sigma);
    tf_pots(1) = tf_pots0;

    res0 = 0;
    D = sigma * z / (A * sigma * (kappa + sigma));

    nu_inf = L * sqrt((alpha/kappa + alpha/sigma)/(Rse));
    D = -L*(kappa*(cosh(nu_inf) - cosh((z-1).*nu_inf)) +...
                         sigma*(1-cosh(z.*nu_inf)+z.*nu_inf.*sinh(nu_inf)))...
                         ./(A*sigma*(kappa + sigma)*...
                         nu_inf.*sinh(nu_inf));

    if electrode == 'pos'
        tf_pots = -tf_pots;
        D = -D;
    end

    if any(isnan(tf_pots))
        error("NAN in tf_pots");
    end
end
