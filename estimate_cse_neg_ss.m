% TODO: Tune ss estimation parameters.
function [A_estimate, B_estimate, C_estimate, D_estimate] = estimate_cse_neg_ss(diffusivity_neg, radius_neg)
    sampling_freq = 190;
    T = 1 / sampling_freq;
    min_T_len = 400;
    num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
    freq_vector = 1 : num_samples - 1;
    s = (2j / T) * tan(pi * freq_vector / num_samples);
    
    num = s / diffusivity_neg * radius_neg ^ 2.0 + 3.0 - 3.0 * radius_neg * sqrt(s / diffusivity_neg) .* coth(radius_neg * sqrt(s / diffusivity_neg));
    den = s .* radius_neg .* (1.0 - radius_neg .* ((s / diffusivity_neg) .^ 0.5) .* coth(radius_neg .* (s / radius_neg) .^ 0.5));
    H_d = num ./ den; % TODO: Add J(z, s)
    h_d = real(ifft(H_d)) * sampling_freq;
    hstep = T * cumsum(h_d);

    step_size = 0.1;
    time_d = 1 : step_size : min_T_len;
    td = T * (1 : num_samples - 1);
    pulse = [0 diff(interp1(td, hstep, time_d))];

    hankel_size = 2;
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

    A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability)];
    B_estimate = [extended_controllability(:, 1)];
    C_estimate = [extended_observability(1, :)];
    D_estimate = [pulse(1)];
end
