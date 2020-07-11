% TODO: Add positive electrode and remove pole at 0.
function [tf_pots, res0, D] = tf_pots(cse, z_coordinates, T_len, sampling_f, electrode, const)
    T = 1 / sampling_f;
    num_samples = 2 ^ (ceil(log2(sampling_f * T_len)));
    f_vector = 1 : num_samples - 0;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    % Declare params
    if electrode == 'neg'
        R_ct = const.R_ct_neg;
        R_film = const.R_film_neg;
        Rse = R_ct + R_film;
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
        as = 3 * eps / Rs;
        Uocv_d = calculate_ocv_derivative_neg(cse, const);
    elseif electrode == 'pos'
        R_ct = const.R_ct_pos;
        R_film = const.R_film_pos;
        Rse = R_ct + R_film;
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
        as = 3 * eps / Rs;
        Uocv_d = calculate_ocv_derivative_pos(cse, const);
    else
        error("Bad electrode selection");
    end

    % Calculate nu
    % z = 0;
    % nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uovc_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));
    % res0 = 0;
    % D = 0;
    % tf_pots = -L * (kappa * (cosh(nu) - cosh(nu * (z - 1))) + sigma * (1 - cosh(z * nu) + z * nu .* sinh(nu))) ./ (A * sigma * (kappa + sigma) * sinh(nu));
    % tf_pots0 = 0;

    z = 0;
    beta = Rs*sqrt(s./Ds);
    cse_j = Rs/Ds./(1-beta.*coth(beta));
    nu = L * sqrt((as/sigma + as/kappa)./(Rse + Uocv_d.*cse_j/F));
    nu_inf = L * sqrt((as/kappa + as/sigma)/(R_ct + R_film));
    tf_pots0 = L*(z-2).*z/(2*A*sigma);
    tf_pots = -L*(kappa*(cosh(nu) - cosh((z-1).*nu)) +...
                      sigma*(1-cosh(z.*nu) + z.*nu.*sinh(nu))) ...
                      ./(A*sigma*(kappa + sigma)*...
                      nu.*sinh(nu));
    D = 10000000000;% -L*(kappa*(cosh(nu_inf) - cosh((z-1).*nu_inf)) +sigma*(1-cosh(z.*nu_inf)+z.*nu_inf.*sinh(nu_inf)))./(A*sigma*(kappa + sigma) * nu_inf.*sinh(nu_inf));
    

    res0 = 0;
    for (i = 1 : size(tf_pots, 2))
        if isnan(tf_pots(1, i))
            tf_pots(1, i) = tf_pots0;
        end
    end

    if electrode == 'pos'
        tf_pots = -tf_pots;
        D = -D;
    end

    if any(isnan(nu_inf))
        error("NAN in nu_inf");
    end
    if any(isnan(D))
        error("NAN in D");
    end
    if any(isnan(tf_pots))
        error("NAN in tf_pots");
    end
end
