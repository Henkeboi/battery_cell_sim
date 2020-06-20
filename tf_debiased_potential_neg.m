% Discretize params
sampling_freq = 400;
T = 1 / sampling_freq;
min_T_len = 100;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
freq_vector = 0 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);
run('calculate_nu_neg.m');

i_app = 0;
TF_pot_neg = zeros(1, size(s, 2));

segments = 1;
TF_pot_neg_vector = zeros(size(s, 2), segments);
for segment = 1 : segments
    for i = 1 : size(s, 2)
        z = segment / segments;
        num = i_app * k_eff_neg * (cosh(nu_neg(i)) -cosh((z - 1) * nu_neg(i))) + o_eff_neg * (1 - cosh(z * nu_neg(i)) + z * nu_neg(i) * sinh(nu_neg(i)));
        den = A_neg * (k_eff_neg + o_eff_neg) * nu_neg(i) * sinh(nu_neg(i));
        z = segments / 100;
        pot_neg_vector(i, segment) = -L_neg * (num / den);
    end
end

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
time_vector = 0 : T_shifted : min_T_len;
A_estimates = [];
B_estimates = [];
C_estimates = [];
D_estimates = [];
for segment = 1 : segments
    segment_pulse = [0 diff(interp1(td, T * cumsum(real(ifft(pot_neg_vector(:, segment))) * sampling_freq), time_vector))];
    hankel_size = 4;
    hankel_matrix = hankel(segment_pulse(:, 2 : end));
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

    A_estimates = [A_estimates; [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)]];
    B_estimates = [B_estimates; [extended_controllability(:, 1); T_shifted]];
    C_estimates = [C_estimates; [extended_observability(1, :), 0]];
    D_estimates = [D_estimates; [segment_pulse(1)]];
    stem(time_vector, segment_pulse, 'filled');
end
