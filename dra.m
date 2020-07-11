function [A_est, B_est, C_est, D_est, T_shifted] = dra(transfer_function, res0, D, sampling_freq, T_len, const)
    % sampling_freq and T_len needs to be the same as in the transfer function.
    T = 1 / sampling_freq;
    num_samples = 2 ^ (ceil(log2(T_len * sampling_freq)));

    tf_c  = real(ifft(transfer_function)) * sampling_freq; 
    h_step = T * cumsum(tf_c); 
    td = T * (0 : num_samples - 1);

    T_shifted = .5;
    time_vector = 0 : T_shifted : T_len;

    h_pulse = [0 diff(interp1(td, h_step, time_vector))];

    hankel_size = 5;
    hankel_matrix = hankel(h_pulse(2 : end));
    H = hankel_matrix(1 : hankel_size, 1 : hankel_size);
    H_shifted = hankel_matrix(2 : hankel_size + 1, 1 : hankel_size);

    [U, S, V] = svds(H);
    system_order = rank(H);
    S = S(1:system_order, 1:system_order);
    U = U(:, 1:system_order);
    V = V(:, 1:system_order);

    sigma = S .^ 0.5;
    extended_observability = U * sigma;
    extended_controllability = sigma * V';

    A_est = [pinv(extended_observability) * H_shifted * pinv(extended_controllability)];
    B_est = [extended_controllability(:, 1)]; 
    C_est = [extended_observability(1, :)];
    D_est = [D];
    res0 = res0(1); % z(1) choosen. TODO: Multiple z or x locations.
end
