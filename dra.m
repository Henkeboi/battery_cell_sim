clearvars;

H_numerator = [1 20 80];
H_denominator = [1 2 8];

sampling_freq = 256;
T = 1 / sampling_freq;
min_T_len = 6.5;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
hd = finely_sampled_continuous_time_impuls_response(H_numerator, H_denominator, sampling_freq, T, num_samples);

h_step = continuous_time_step_response(hd, T);

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
[time_vector, h_pulse] = unit_pulse_response(h_step, T_shifted, min_T_len, td);
[A_est, B_est, C_est, D_est] = ho_kalman(h_pulse, 2)

% Plot
stem(time_vector, h_pulse,'filled'); 
hold on;
H_true = tf([1 20 80], [1 2 8]);
[himpDiscTrue, timpDiscTrue] = impulse(c2d(H_true, T_shifted), 5);
himpDiscTrue = T_shifted * himpDiscTrue;
plot(timpDiscTrue, himpDiscTrue, 'r.', 'markersize',8);
axis([-0.01 5 -0.8 2.3]);

function [time_vector, h_pulse] = unit_pulse_response(h_step, T_shifted, min_T_len, td)
    time_vector = 0 : T_shifted : min_T_len;
    h_pulse = [0 diff(interp1(td, h_step, time_vector))];
end

function [h_step] = continuous_time_step_response(hd, T)
    h_step = T * cumsum(hd); 
end

function [hd] = finely_sampled_continuous_time_impuls_response(H_numerator, H_denominator, sampling_freq, T, num_samples)
    freq_vector = 0 : num_samples - 1;
    s = (2j / T) * tan(pi * freq_vector / num_samples);
    % Hd[f] is the DFT of H(z) found by the bilinear transform of H(s) which is the continuous time impulse response
    Hd = (s .^2 + 20 * s + 80) ./ (s .^ 2 + 2 * s + 8); 
    % Approximation to the continuous time impulse response
    hd = real(ifft(Hd)) * sampling_freq; 
end

function [A_estimate, B_estimate, C_estimate, D_estimate] = ho_kalman(unit_pulse, hankel_size)
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

    A_estimate = pinv(extended_observability) * H_shifted * pinv(extended_controllability);
    B_estimate = extended_controllability(:, 1);
    C_estimate = extended_observability(1, :);
    D_estimate = unit_pulse(1);
end
