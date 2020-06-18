% Discretize params 
sampling_freq = 2560;
T = 1 / sampling_freq;
min_T_len = 10;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
freq_vector = 0 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);
run('calculate_nu_neg.m');

z = 1; % Random z-coordinate
TF_pot_neg = zeros(1, size(s, 2));
for i = 1 : size(s, 2)
    num = k_eff_neg * (cosh(nu_neg(i)) -cosh((z - 1) * nu_neg(i))) + o_eff_neg * (1 - cosh(z * nu_neg(i)) + z * nu_neg(i) * sinh(nu_neg(i)));
    den = A_neg * (k_eff_neg + o_eff_neg) * nu_neg(i) * sinh(nu_neg(i));
    TF_pot_neg(i) = -L_neg * (num / den);
end

tf_pot_neg = real(ifft(TF_pot_neg)) * sampling_freq;
pot_neg_step = T * cumsum(tf_pot_neg);

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
time_vector = 0 : T_shifted : min_T_len;
pot_neg_pulse = [0 diff(interp1(td, pot_neg_step, time_vector))];

hankel_size = 2;
hankel_matrix = hankel(pot_neg_pulse(2:end));
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
C_est = [extended_observability(1, :), 0];
D_est = [pot_neg_pulse(1)];
disp(A_est)

% Plot
stem(time_vector, pot_neg_pulse,'filled');
grid on;
axis([-0.01 5 -0.8 2.3]);
