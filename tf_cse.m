function [tf_cse, res0, D] = tf_cse(cse, z_coordinates, T_len, sampling_f, electrode, const)
    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);

    z = 0.1;

    tf_cse = (sigma * cosh(nu .* z) + kappa * cosh(nu * (z - 1))) ./ (as * F * L * A * Ds * (kappa + sigma) .* sinh(nu)) ...
        .* (Rs .* nu ./ (1 - beta .* coth(beta)));
   
    tf_cse(1, 1) = 0;
    for i = 1 : size(tf_cse, 1)
        for j = 1 : size(tf_cse, 2)
            if isnan(tf_cse(i, j))
                if z == 0
                    tf_cse(i, j)  = kappa ./ (as * F * L * A * Ds * (kappa + sigma)) .* (Rs * nu(i, j) / (1 - beta(i, j) * coth(beta(i, j))));
                elseif z == 1
                    tf_cse(i, j) = sigma ./ (as * F * L * A * Ds * (kappa + sigma)) .* (Rs * nu(i, j) / (1 - beta(i, j) * coth(beta(i, j))));
                else
                    tf_cse(i, j) = 0;
                end
            end
        end
    end

    res0 = 1 ./ (eps * A * F * L .* s);
    D = 0;

    if electrode == 'pos'
        tf_cse = -tf_cse;
        res0 = -res0;
    end

    if any(isnan(tf_cse))
        error("NAN In tf_cse");
    end

end
