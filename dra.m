clearvars;

% G_numerator = [1 20 80];
% G_denominator = [1 2 8];

sampling_freq = 256;
T = 1 / sampling_freq;
min_T_len = 6.5;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
hd = impuls_response(sampling_freq, T, num_samples);
% hd = finely_sampled_continuous_time_impuls_response(G_numerator, G_denominator, sampling_freq, T, num_samples);
h_step = continuous_time_step_response(hd, T);

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
[time_vector, h_pulse] = unit_pulse_response(h_step, T_shifted, min_T_len, td);
[A_est, B_est, C_est, D_est] = ho_kalman(h_pulse, 4, T_shifted, 0);

% Plot
stem(time_vector, h_pulse,'filled');
hold on;
H_true = tf([-0.5 -3], [1 6 2]);
[himpDiscTrue, timpDiscTrue] = impulse(c2d(H_true, T_shifted), 5);
himpDiscTrue = T_shifted * himpDiscTrue;
plot(timpDiscTrue, himpDiscTrue, 'r.', 'markersize', 8);
grid on;
axis([-0.01 5 -0.8 2.3]);

function [time_vector, h_pulse] = unit_pulse_response(h_step, T_shifted, min_T_len, td)
    time_vector = 0 : T_shifted : min_T_len;
    h_pulse = [0 diff(interp1(td, h_step, time_vector))];
end

function [h_step] = continuous_time_step_response(hd, T)
    h_step = T * cumsum(hd); 
end

function [hd] = impuls_response(sampling_freq, T, num_samples)
    freq_vector = 0 : num_samples - 1;
    s = (2j / T) * tan(pi * freq_vector / num_samples);
    Hd = (-0.5 * s - 3) ./ (s .^ 2 + 6 * s + 2);
    hd = real(ifft(Hd)) * sampling_freq;
end

function [hd] = finely_sampled_continuous_time_impuls_response(G_numerator, G_denominator, sampling_freq, T, num_samples)
    freq_vector = 0 : num_samples - 1;
    s = (2j / T) * tan(pi * freq_vector / num_samples);

    G_numerator_accumulated = 0;
    for i = 1 : size(G_numerator, 2)
        G_numerator_accumulated = G_numerator_accumulated + s .^ (size(G_numerator, 2) - i) * G_numerator(i);
    end
    G_denominator_accumulated = 0;
    for i = 1 : size(G_denominator, 2)
        G_denominator_accumulated = G_denominator_accumulated + s .^ (size(G_denominator, 2) - i) * G_denominator(i);
    end

    % Hd[f] is the DFT of H(z) found by the bilinear transform of H(s) which is the continuous time impulse response.
    Hd = G_numerator_accumulated ./ G_denominator_accumulated; 
    % Approximation to the continuous time impulse response.
    hd = real(ifft(Hd)) * sampling_freq; 
end

function [A_estimate, B_estimate, C_estimate, D_estimate] = ho_kalman(unit_pulse, hankel_size, T_shifted, res0)
    hankel_matrix = hankel(unit_pulse(2:end));
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

    A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)];
    B_estimate = [extended_controllability(:, 1); T_shifted]; 
    C_estimate = [extended_observability(1, :), res0];
    D_estimate = [unit_pulse(1)];
end
