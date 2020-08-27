function [tf_potse, res0, D] = tf_j(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    z = 0.5;
    tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    tf_potse0 = 10e10;%L * (z * z * (kappa + sigma) / 2 + kappa * (1 / 2 - z));
    tf_potse(1) = tf_potse0;

    res0 = 0;
    D = 0;

    if electrode == 'pos'
        tf_potse = -tf_potse;
    end
   
    if any(isnan(tf_potse))
        error("NAN in tf_potse");
    end
end
