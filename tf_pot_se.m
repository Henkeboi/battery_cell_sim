function [tf_phi] = j_tf(const)
    sampling_f = 200;
    T = 1 / sampling_f;
    min_T_len = 5;
    num_samples = 2 ^ (ceil(log2(sampling_f * min_T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end

    % Params
    cse_neg = 1;
    z = 0.5; % z coordinate

    Uovc_d = calculate_ocv_derivative_neg(cse_neg, const);

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
    nu_num = L * sqrt((alpha / sigma) + (alpha / kappa));
    nu_den = sqrt(Rse + (Rs * Uovc_d) / (F * D)) * sqrt(1 ./ (1 - beta .* coth(beta)));
    nu = nu_num ./ nu_den;
    tf_phi = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu)); % PHI / Iapp 
    tf_phi0 = (L * z * z * (kappa + sigma) - 2 * L * z * kappa + L * kappa) / (2 * A * kappa * sigma); % Solved with Maple when s -> 0
    nan_indexes = isnan(tf_phi);
    tf_phi(nan_indexes) = tf_phi0; % Replace with tf_phi0.

    hd = real(ifft(tf_phi)) * sampling_f;
    disp(hd)
    

end
