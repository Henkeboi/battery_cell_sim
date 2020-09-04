function [tf_delta_cse, res0, D] = tf_delta_cse(cse, z, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    tf_delta_cse = nu .* (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ...
                    ./ (as * F * L * A * (kappa + sigma) * sinh(nu)) ...
                    .* (s .* Rs ^ 2 / Ds + 3 - 3 * Rs .* sqrt(s ./ Ds) .* coth(Rs .* sqrt(s ./ Ds))) ...
                    ./ (s .* Rs .* (1 - Rs * sqrt(s / Ds) .* coth(Rs .* sqrt(s ./ Ds))));


    tf_delta_cse(1, 1) = -1 / (5 * Rs) * 1 / (as + F * L * A);
    res0 = 0;
    D = NaN;
    if electrode == 'pos'
        tf_delta_cse = -tf_delta_cse;
        res0 = -res0;
    end

    for i = 1 : size(tf_delta_cse, 1)
        for j = 1 : size(tf_delta_cse, 2)
            if isnan(tf_delta_cse(i, j))
                nu = nu(i, j);
                s = s(i, j);
                disp(nu .* (sigma * cosh(nu * z) + kappa * cosh(nu * (z - 1))) ...
                    ./ (as * F * L * A * (kappa + sigma) * sinh(nu)) ...
                    .* (s .* Rs ^ 2 / Ds + 3 - 3 * Rs .* sqrt(s ./ Ds) .* coth(Rs .* sqrt(s ./ Ds))) ...
                    ./ (s .* Rs .* (1 - Rs * sqrt(s / Ds) .* coth(Rs .* sqrt(s ./ Ds)))));
                disp(i + " " + j)
            end
        end
    end

    if any(isnan(tf_delta_cse))
        error("NAN In tf_delta_cse");
    end
end
