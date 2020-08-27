function [tf_cse, res0, D] = tf_cse(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    z = 0.0000001; % z = 0 not good with current impl.
    nu = L * sqrt(as / sigma + as / kappa) ./ sqrt((Rse + Rs * Uocv_d / (F * Ds)) ./ (1 - beta .* coth(beta)));
    %[M, index] = max(nu);
    %nu_inf = L * sqrt((as / sigma + as / kappa) / Rse);
    %nu(index) = nu_inf;

    %tf_potse = L * (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ./ (A * sigma * kappa * nu .* sinh(nu));
    %tf_potse0 = L * (z * z * (kappa + sigma) / 2 + kappa * (1 / 2 - z));
    %tf_potse(1) = tf_potse0;
    %tf_j = nu .* nu .* tf_potse / (as * F * L * L * (1 / kappa + 1 / sigma));
    %tf_j(1) = 0.0;
    %tf_cse = Rs / Ds ./ (1 - beta .* coth(beta)) .* tf_j;
    %tf_cse(1) = Rs / Ds * tf_j(1);
    tf_cse = Rs * sigma * nu .* coth(nu * z) ./ (as * F * L * A * Ds * (kappa + sigma) * (1 - beta .* coth(beta)));
    tf_cse = tf_cse + Rs * kappa * nu .* (exp(z) ./ (exp(nu) - 1) + exp(-z * nu) ./ (1 - exp(-nu))) ./ (as * F * L * A * Ds * (kappa + sigma) * (1 - beta .* coth(beta)));
    tf_cse(1) = 0;

    res0 = 0;
    D = 0;

    if electrode == 'pos'
        tf_cse = -tf_cse;
        res0 = -res0;
    end

    if any(isnan(tf_cse))
        error("NAN In tf_cse");
    end
end
