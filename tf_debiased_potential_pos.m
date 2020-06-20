% Discretize params
sampling_freq = 150;
T = 1 / sampling_freq;
min_T_len = 1;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
freq_vector = 0 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);
run('calculate_nu_pos.m');

i_app = 1;
TF_pot_pos = zeros(1, size(s, 2));

segments = 1;
TF_pot_pos_vector = zeros(size(s, 2), segments);
for segment = 1 : segments
    for i = 1 : size(s, 2)
        z = (segment - 1) / segments;
        % NAN error from Matlabs cosh.
        cosh_on_z_1_nu_pos = (exp(nu_pos(i)) + exp(-nu_pos(i)) -exp((z - 1) * nu_pos(i)) + exp(-(z - 1) * nu_pos(i))) / 2;
        cosh_on_z_nu_pos = (exp(z * nu_pos(i)) + exp(z * nu_pos(i))) / 2;
        sinh_on_nu_pos = (exp(nu_pos(i)) - exp(-nu_pos(i))) / 2;
        disp(sinh_on_nu_pos);
        num = i_app * (k_eff_pos * cosh_on_z_1_nu_pos + o_eff_pos * (1 - cosh_on_z_nu_pos + (z * nu_pos(i) * sinh_on_nu_pos)));
        den = A_pos * (k_eff_pos + o_eff_pos) * nu_pos(i) * sinh(nu_pos(i));
        % if (num / den ~= num / den)
        %     disp("a");
        % end
        TF_pot_pos_vector(i, segment) = L_neg * (num / den);
    end
end
% for i = 1 : size(pot_pos_vector, 1)
%     if pot_pos_vector(i, 1) ~= pot_pos_vector(i, 1)
%         disp(i)
%         disp(pot_pos_vector(:, 1));
%     end
% end

td = T * (0 : num_samples - 1);
T_shifted = 0.1;
time_vector = 0 : T_shifted : min_T_len;
A_estimates = [];
B_estimates = [];
C_estimates = [];
D_estimates = [];
for segment = 1 : segments
    segment_pulse = [0 diff(interp1(td, T * cumsum(real(ifft(TF_pot_pos_vector(:, segment))) * sampling_freq), time_vector))];
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
