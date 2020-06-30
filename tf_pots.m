% TODO: Add positive electrode and remove pole at 0.
function [tf_pots, res0, D, sampling_f, T_len] = tf_pots(cse, z_coordinates, const, electrode)
    sampling_f = 200;
    T = 1 / sampling_f;
    T_len = 4;
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

    % Calculate nu
    z = 0.01;
    nu = L * sqrt(alpha / sigma + alpha / kappa) ./ sqrt(Rse + Uovc_d / F / Ds * (1 ./ (Rs * beta .* coth(beta))));
    res0 = 0;
    D = 0;

    % Calculate TF potential in solid.
    tf_pots = -L * (kappa * (cosh(nu) - cosh(nu * (z - 1))) + sigma * (1 - cosh(z * nu) + z * nu .* sinh(nu))) ./ (A * sigma * (kappa + sigma) * sinh(nu));
    tf_pots0 = 0; % TODO: Find this value.
    for (i = 1 : size(tf_pots, 2))
        if isnan(tf_pots(1, i)) && s(1, i) == 0
            tf_pots(1, i) = tf_pots0;
        end
    end

    if electrode == 'pos'
        tf_pots = -tf_pots;
    end

    if any(isnan(tf_pots))
        error("NAN in tf_j");
    end
end
