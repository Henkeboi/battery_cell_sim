function [s, nu] = calculate_nu(cse, T_len, sampling_f, electrode, const)
    T = 1 / sampling_f;
    num_samples = 2 ^ (ceil(log2(sampling_f * T_len)));
    f_vector = 0 : num_samples - 1;
    s = zeros(1, size(f_vector, 2));
    for i = 1 : size(f_vector, 2)
        s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    end
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d] = get_electrode_constants(s, electrode, const);

    nu = L * sqrt((as / sigma + as / kappa) ./ (Rse + Rs * Uocv_d ./ (F * Ds - F * Ds * beta .* coth(beta))));
    nu(1) = 0.0;

    if any(isnan(nu))
        error("NAN in calculate nu")
    end
end
