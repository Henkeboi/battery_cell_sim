function [j_neg_ss] = estimate_li_flux_ss(cse_neg, const)
    F = const.F;
    L_neg = const.L_neg;
    A_neg = const.A_neg;
    resistivity_neg = const.resistivity_neg;
    solid_max_c_neg = const.solid_max_c_neg;
    k_eff_neg = const.k_eff_neg;
    sigma_eff_neg = const.sigma_eff_neg;
    diffusivity_neg = const.diffusivity_neg;
    asymmetric_charge_transfer_neg = const.asym_ct_neg;
    R_neg = const.R_neg;
    R_solid_electrolyte_neg = const.R_solid_electrolyte_neg;
    ionic_conductivity_neg = const.ionic_conductivity_neg;

    sampling_freq = 300;
    T = 1 / sampling_freq;
    min_T_len = 10;
    num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
    freq_vector = 0 : num_samples;
    s = zeros(1, size(freq_vector, 2));
    for i = 1 : size(freq_vector, 2)
        s(i) = 2j * tan(pi * freq_vector(i));
    end

    ocv_derivative_neg = calculate_ocv_derivative_neg(cse_neg, const);
    nu_neg = calculate_nu_neg(s, ocv_derivative_neg, const);

    H_d = zeros(1, size(s, 2));
    res0 = 1 / (L_neg * asymmetric_charge_transfer_neg * F * A_neg);
    for i = 1 : size(s, 2)
        z_coordinate = 0.5;
        numerator = L_neg * sqrt(asymmetric_charge_transfer_neg / (resistivity_neg + ionic_conductivity_neg));
        denominator_1 = R_solid_electrolyte_neg;
        denominator_2 = R_neg * ocv_derivative_neg / (F * diffusivity_neg);
        denominator_3 = 1.0;
        denominator_4 = R_solid_electrolyte_neg * sqrt(s(i) / diffusivity_neg);
        denominator_5 = (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) + 1.0) / (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) - 1.0);
        denominator = sqrt(denominator_1 + denominator_2 * (1.0 / (denominator_3 - denominator_4 * denominator_5)));
        num = nu_neg(i) * (sigma_eff_neg * cosh(nu_neg(i) .* z_coordinate) + k_eff_neg * cosh(nu_neg(i) * (z_coordinate - 1)));
        den = asymmetric_charge_transfer_neg .* F .* L_neg .* A_neg .* (k_eff_neg + sigma_eff_neg) .* sinh(nu_neg(i));
        H_d(i) = num / den;
    end
    H_d = H_d - res0;
    omega = logspace(-1, 3, size(H_d, 2)); % create freq. axis in rad/sec
    s = 1j * omega;

    h_d = real(ifft(H_d)) * sampling_freq;
    t_d = T * (1 : size(H_d, 2));

    tran = tf(H_d);

    step_size = 0.5;
    time_d = 1 : step_size : min_T_len;
    td = T * (1 : num_samples - 1);
    h_step = T * cumsum(h_d);
    pulse = [0 diff(interp1(td, h_step, time_d))];

    hankel_size = 3;
    hankel_matrix = hankel(pulse(:, 2 : end));
    H = hankel_matrix(1:hankel_size, 1:hankel_size);
    H_shifted = hankel_matrix(2:hankel_size + 1, 1:hankel_size);
    [U, S, V] = svds(H);
    system_order = rank(H);
    S = S(1:system_order, 1:system_order);
    U = U(:, 1:system_order);
    V = V(:, 1:system_order);
    sigma = S .^ 0.5;
    extended_observability = U * sigma;
    extended_controllability = sigma * V';

    %A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)];
    %B_estimate = [extended_controllability(:, 1); 0.1];
    %C_estimate = [extended_observability(1, :), 0];
    %D_estimate = [pulse(1)];

    A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability)];
    B_estimate = [extended_controllability(:, 1)];
    C_estimate = [extended_observability(1, :)];
    D_estimate = [pulse(1)];
    disp(eig(A_estimate))

    j_neg_ss = ss(A_estimate, B_estimate, C_estimate, D_estimate, step_size);

    % eig_values = eig(A_estimate);
    % P = [-1 eig_values(2) eig_values(3)];
    % K = place(A_estimate, B_estimate, P);
    % Acl = A_estimate - B_estimate * K;
    % sysctl = ss(Acl, B_estimate, C_estimate, D_estimate);
    % Kdc = dcgain(sysctl);
    % Kr = 1 / Kdc;
    % sysctl_scaled = ss(Acl, Kr * B_li_neg, C_li_neg, D_li_neg);
    % A_li_neg = Acl;
    % disp(eig(Acl))
    % B_li_neg = Kr * B_li_neg;
    % step(sysctl_scaled)
    % impulse(sysctl_scaled);
    % bode(sysctl_scaled);
end
