function [A_estimate, B_estimate, C_estimate, D_estimate] = estimate_li_flux_ss(cse_neg, solid_max_concentration_neg, sigma_eff_neg, k_eff_neg, L_neg, asymmetric_charge_transfer_neg, resistivity_neg, ionic_conductivity_neg, R_solid_electrolyte_neg, R_neg, diffusivity_neg, F, A_neg)
    sampling_freq = 200000;
    T = 1 / sampling_freq;
    min_T_len = 4;
    num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
    freq_vector = 1 : num_samples - 1;
    s = (2j / T) * tan(pi * freq_vector / num_samples);

    ocv_derivative_neg = calculate_ocv_derivative_neg(cse_neg, solid_max_concentration_neg);
    %nu_neg = calculate_nu_neg(s, ocv_derivative_neg, sigma_eff_neg, k_eff_neg, L_neg, asymmetric_charge_transfer_neg, resistivity_neg, ionic_conductivity_neg, R_solid_electrolyte_neg, R_neg, diffusivity_neg, F);
     
    charge_cycle = zeros(1, size(s, 2));
    for i = 1 : size(s, 2)
        charge_cycle(i) = i;
    end
    H_d = zeros(1, size(s, 2));
    for i = 1 : size(s, 2)
        z_coordinate = 0.5;
        numerator = L_neg * sqrt(asymmetric_charge_transfer_neg / (resistivity_neg + ionic_conductivity_neg));
        denominator_1 = R_solid_electrolyte_neg;
        denominator_2 = R_neg * ocv_derivative_neg / (F * diffusivity_neg);
        denominator_3 = 1.0;
        denominator_4 = R_solid_electrolyte_neg * sqrt(s(i) / diffusivity_neg);
        denominator_5 = (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) + 1.0) / (exp(2.0 * R_neg / sqrt(s(i) / diffusivity_neg)) - 1.0);
        denominator = sqrt(denominator_1 + denominator_2 * (1.0 / (denominator_3 - denominator_4 * denominator_5)));
        nu_neg = numerator / denominator;
        num = nu_neg * (sigma_eff_neg * cosh(nu_neg .* z_coordinate) + k_eff_neg * cosh(nu_neg * (z_coordinate - 1)));
        den = asymmetric_charge_transfer_neg .* F .* L_neg .* A_neg .* (k_eff_neg + sigma_eff_neg) .* sinh(nu_neg);
        H_d(1, i) = charge_cycle(i) * (num) ./ (den);
    end

    z_coordinate = 0.5;
    num = nu_neg .* (sigma_eff_neg .* cosh(nu_neg .* z_coordinate) + k_eff_neg .* cosh(nu_neg .* (z_coordinate - 1)));
    den = asymmetric_charge_transfer_neg .* F .* L_neg .* A_neg .* (k_eff_neg + sigma_eff_neg) .* sinh(nu_neg);
    
    omega = logspace(-1,3,size(s, 2)); % create freq. axis in rad/sec
    s = 1j*omega; % create s = j*omega
    H = (s.^2+20*s+80)./(s.^2+2*s+8); % compute cplx. freq. response
    semilogx(omega,20*log10(abs(H_d))); %

    h_d = real(ifft(H_d)) * sampling_freq;
    td = T * (1 : size(H_d, 2));

    % plot(td, h_d, 'bx', 'markersize', 1);
    
    step_size = 0.1;
    time_d = 1 : step_size : min_T_len;
    td = T * (1 : num_samples - 1);
    h_step = T * cumsum(h_d);
    pulse = [0 diff(interp1(td, h_step, time_d))];


    A_estimate = [];
    B_estimate = [];
    C_estimate = [];
    D_estimate = [];


    % hankel_size = 2;
    % hankel_matrix = hankel(pulse(:, 2 : end));
    % H = hankel_matrix(1:hankel_size, 1:hankel_size);
    % H_shifted = hankel_matrix(2:hankel_size + 1, 1:hankel_size);
    % [U, S, V] = svds(H);
    % system_order = rank(H);
    % S = S(1:system_order, 1:system_order);
    % U = U(:, 1:system_order);
    % V = V(:, 1:system_order);
    % sigma = S .^ 0.5;
    % extended_observability = U * sigma;
    % extended_controllability = sigma * V';

    % A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability)];
    % B_estimate = [extended_controllability(:, 1)];
    % C_estimate = [extended_observability(1, :)];
    % D_estimate = [pulse(1)];
end
