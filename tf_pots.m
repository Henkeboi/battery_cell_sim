function [tf_pots, res0, D] = tf_pots(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    z = 0.0;
    tf_pots = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    tf_pots0 = -L * z * (1 - z / 2) / (A * sigma);
    tf_pots(1) = tf_pots0;

    res0 = 0;
    D = sigma * z / (A * sigma * (kappa + sigma));

    if electrode == 'pos'
        tf_pots = -tf_pots;
        D = -D;
    end

    if any(isnan(tf_pots))
        disp(find(isnan(tf_pots)))
        error("NAN in tf_pots");
    end
end
