function [tf_j, res0, D] = tf_j(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    z = 0.4;
    tf_j = nu .* (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (as * A * F * L * (kappa + sigma) * sinh(nu));
    tf_j(1) = 1 / (as * F * L * A);
    [row, col] = find(isnan(tf_j));
    tf_j(row, col) = 0;
    res0 = 0;
    D = NaN;

    if electrode == 'pos'
        tf_j = -tf_j;
    end
   
    if any(isnan(tf_j))
        error("NAN in tf_j");
    end
end
