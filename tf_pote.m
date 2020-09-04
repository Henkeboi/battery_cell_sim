function [tf_pote, res0, D] = tf_pote(cse, x, T_len, sampling_f, electrode, const)
    %TODO:
    ce = const.ce0_neg;
    % END

    if x <= const.L_neg
        region = 'neg';
        electrode = 'neg';
        [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
        [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);
    elseif x < const.L_neg + const.L_sep
        region = 'sep';
        electrode = 'neg';
        L_sep = const.L_sep;
        conductivity = 4.1253 * 10e-2 + 5.007 * 10e-4 * ce - 4.7212 * 10e-7 * ce ^ 2 + 1.5094 * 10e-10 * ce ^ 3 - 1.6018 * 10e-14 * ce ^ 4;
        kappa_sep = conductivity * const.eps_e_sep ^ const.brug_sep;
        [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
        [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);
    elseif x == const.L_neg + const.L_sep
        region = 'bnd';
        electrode = 'neg';
        L_sep = const.L_sep;
        conductivity = 4.1253 * 10e-2 + 5.007 * 10e-4 * ce - 4.7212 * 10e-7 * ce ^ 2 + 1.5094 * 10e-10 * ce ^ 3 - 1.6018 * 10e-14 * ce ^ 4;
        kappa_sep = conductivity * const.eps_e_sep ^ const.brug_sep;
        [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
        [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);
    elseif x < const.L_neg + const.L_sep + const.L_pos
        region = 'pos';
        L_sep = const.L_sep;
        conductivity = 4.1253 * 10e-2 + 5.007 * 10e-4 * ce - 4.7212 * 10e-7 * ce ^ 2 + 1.5094 * 10e-10 * ce ^ 3 - 1.6018 * 10e-14 * ce ^ 4;
        kappa_sep = conductivity * const.eps_e_sep ^ const.brug_sep;
        [s_n, nu_n] = calculate_s_nu(cse, T_len, sampling_f, 'neg', const);
        [Rct_n, Rfilm_n, Rse_n, Rs_n, as_n, L_n, F_n, Ds_n, A, alpha_n, sigma_n, kappa_n, beta_n, eps_n, cs0_n, Uocv_d_n, transfer_number_n] = get_electrode_constants(s_n, 'neg', const);
        [s_p, nu_p] = calculate_s_nu(cse, T_len, sampling_f, 'pos', const);
        [Rct_p, Rfilm_p, Rse_p, Rs_p, as_p, L_p, F_p, Ds_p, A, alpha_p, sigma_p, kappa_p, beta_p, eps_p, cs0_p, Uocv_d_p, transfer_number_p] = get_electrode_constants(s_p, 'pos', const);
        L_tot = L_n + L_sep + L_n;
    elseif x == const.L_neg + const.L_sep + const.L_pos
        region = 'end';
        L_sep = const.L_sep;
        conductivity = 4.1253 * 10e-2 + 5.007 * 10e-4 * ce - 4.7212 * 10e-7 * ce ^ 2 + 1.5094 * 10e-10 * ce ^ 3 - 1.6018 * 10e-14 * ce ^ 4;
        kappa_sep = conductivity * const.eps_e_sep ^ const.brug_sep;
        [s_n, nu_n] = calculate_s_nu(cse, T_len, sampling_f, 'neg', const);
        [Rct_n, Rfilm_n, Rse_n, Rs_n, as_n, L_n, F_n, Ds_n, A, alpha_n, sigma_n, kappa_n, beta_n, eps_n, cs0_n, Uocv_d_n, transfer_number_n] = get_electrode_constants(s_n, 'neg', const);
        [s_p, nu_p] = calculate_s_nu(cse, T_len, sampling_f, 'pos', const);
        [Rct_p, Rfilm_p, Rse_p, Rs_p, as_p, L_p, F_p, Ds_p, A, alpha_p, sigma_p, kappa_p, beta_p, eps_p, cs0_p, Uocv_d_p, transfer_number_p] = get_electrode_constants(s_p, 'pos', const);
        L_tot = L_n + L_sep + L_n;
    else
        error("x location to large.")
    end
    
    if region == 'neg' % At the negative electrode / separator boundary.
        tf_pote = -L * (sigma - kappa) * tanh(nu ./ 2) ./ (A * kappa * (kappa + sigma) .* nu) - L / (A * (kappa + sigma));
        tf_pote(1, 1) = -L * (sigma - kappa) / (2 * A * kappa * (kappa + sigma)) - L / (A * (kappa + sigma));
        res0 = 0;
        D = - L / (A * (kappa + sigma));
    elseif region == 'sep'
        tf_pote = -L * (sigma - kappa) * tanh(nu ./ 2) ./ (A * kappa * (kappa + sigma) .* nu) - L / (A * (kappa + sigma)) - (x - L) / (A * kappa_sep);
        tf_pote(1, 1) = -L * (sigma - kappa) / (2 * A * kappa * (kappa + sigma)) - L / (A * (kappa + sigma)) - (x - L) / (A * kappa_sep);
        res0 = 0;
        D = - L / (A * (kappa + sigma)) - (x - L) / (A * kappa_sep);
    elseif region == 'bnd' % At the separator / positive electrode boundary.
        tf_pote = -L * (sigma - kappa) * tanh(nu ./ 2) ./ (A * kappa * (kappa + sigma) .* nu) - L / (A * (kappa + sigma)) - L_sep / (A * kappa_sep);
        tf_pote(1, 1) = -L * (sigma - kappa) / (2 * A * kappa * (kappa + sigma)) - L / (A * (kappa + sigma)) - L_sep / (A * kappa_sep);
        res0 = 0;
        D = - L / (A * (kappa + sigma)) - L_sep / (A * kappa_sep);
    elseif region == 'pos'
        tf_pote = - L_n * (sigma_n - kappa_n) * tanh(nu_n ./ 2) ./ (A * kappa_n * (kappa_n + sigma_n) .* nu_n) ...
                  - L_n ./ (A * (kappa_n + sigma_n)) ...
                  - L_sep / (A * kappa_sep) ...
                  - L_p * (1 - cosh((L_n + L_sep - x) .* nu_p / L_p)) ./ (A * (kappa_p + sigma_p) .* sinh(nu_p) .* nu_p) ...
                  - L_p * sigma_p * (cosh(nu_p) - cosh((L_tot - x) / L_p .* nu_p)) ./ (A * kappa_p * (kappa_p + sigma_p) .* sinh(nu_p) .* nu_p) ...
                  - (x - L_n - L_sep) ./ (A * (kappa_p + sigma_p));

        tf_pote(1, 1) = -L_n * (sigma_n - kappa_n) ./ (A * kappa_n * (kappa_n + sigma_n)) ...
            - L_n / (A * (kappa_n + sigma_n)) ...
            - L_sep / (A * kappa_sep) ...
            + L_p * (((L_n + L_sep - x) / L_p) ^ 2) / (2 * A * (kappa_p + sigma_p)) ...
            - L_p * sigma_p * (1 - (((L_tot - x) / L_p) ^ 2)) / (2 * A * kappa_p * (kappa_p + sigma_p)) ...
            - (x - L_n - L_sep) / (A * (kappa_p + sigma_p));

        for i = 1 : size(tf_pote, 1)
            for j = 1 : size(tf_pote, 2)
                if isnan(tf_pote(i, j))
                    if isinf(exp(nu_p(j)))
                        tf_pote(i, j) = - L_n * (sigma_n - kappa_n) * tanh(nu_n(j) ./ 2) ./ (A * kappa_n * (kappa_n + sigma_n) .* nu_n(j)) ...
                                        - L_n ./ (A * (kappa_n + sigma_n)) ...
                                        - L_sep / (A * kappa_sep) ...
                                        - L_p * sigma_p * coth(nu_p(j)) ./ (A * kappa_p * (kappa_p + sigma_p) .* nu_p(j)) ...
                                        - (x - L_n - L_sep) ./ (A * (kappa_p + sigma_p));
                    end
                end
            end
        end
        res0 = 0;
        D = -L_n / (A * (kappa_n + sigma_n)) - L_sep / (A * kappa_sep) - L_p * - (x - L_n - L_sep) / (A * (kappa_p + sigma_p));
    elseif region == 'end'
        tf_pote = - L_n * (sigma_n - kappa_n) * tanh(nu_n ./ 2) ./ (A * kappa_n * (kappa_n + sigma_n) .* nu_n)...
                    - L_n / (A * (kappa_n + sigma_n)) ...
                    - L_sep / (A * kappa_sep) ...
                    - L_p * (sigma_p - kappa_p) * tanh(nu_p ./ 2) ./ (A * kappa_p * (kappa_p + sigma_p) .* nu_p) ...
                    - L_p / (A * (kappa_p + sigma_p));
        
        tf_pote(1, 1) = - L_n * (sigma_n - kappa_n) / (2 * A * kappa_n * (kappa_n + sigma_n)) ...
                        - L_n / (A * (kappa_n + sigma_n)) ...
                        - L_sep / (A * kappa_sep) ...
                        - L_p * (sigma_p - kappa_p) / (2 * kappa_p * (kappa_p + sigma_p)) ...
                        - L_p / (A * (kappa_p + sigma_p));
        res0 = 0;
        D = - L_n / (A * (kappa_n + sigma_n)) ...
            - L_sep / (A * kappa_sep) ...
            - L_p / (A * (kappa_p + sigma_p));
    end

    if any(isnan(tf_pote)) 
        error("NaN in tf pote")
    end
end
