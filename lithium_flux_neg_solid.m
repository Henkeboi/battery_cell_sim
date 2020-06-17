clearvars;
run('parameters.m')

% Params
sampling_freq = 256;
T = 1 / sampling_freq;
min_T_len = 6.5;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));

% Calculate hd
freq_vector = 0 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);
Hd = (-0.5 * s - 3) ./ (s .^ 2 + 6 * s + 2);
hd = real(ifft(Hd)) * sampling_freq;

% Calculate h step
h_step = T * cumsum(hd); 

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
time_vector = 0 : T_shifted : min_T_len;
h_pulse = [0 diff(interp1(td, h_step, time_vector))];

freq_vector = 0 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);
run('calculate_nu_neg.m');
disp(nu_neg);

% Hd[f] is the DFT of H(z) found by the bilinear transform of H(s) which is the continuous time impulse response.
Hd = (0.5 * s - 3) ./ (s .^ 2 + 6 * s + 2);
% Approximation to the continuous time impulse response.
hd = real(ifft(Hd)) * sampling_freq; 
res0 = 0;
hankel_size = 2;
hankel_matrix = hankel(h_pulse(2:end));
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

A_est = [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)];
B_est = [extended_controllability(:, 1); T_shifted]; 
C_est = [extended_observability(1, :), res0];
D_est = [h_pulse(1)];

% Plot
stem(time_vector, h_pulse,'filled');
hold on;
H_true = tf([-0.5 -3], [1 6 2]);
[himpDiscTrue, timpDiscTrue] = impulse(c2d(H_true, T_shifted), 5);
himpDiscTrue = T_shifted * himpDiscTrue;
plot(timpDiscTrue, himpDiscTrue, 'r.', 'markersize', 8);
grid on;
axis([-0.01 5 -0.8 2.3]);
