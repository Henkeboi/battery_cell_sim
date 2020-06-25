clearvars;
run('parameters.m');

cse_pos = 10000;
cse_neg = 10000;
z_coordinates = [0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

[all_tf_j_pos, res0, D, sampling_freq, T_len] = tf_j(cse_pos, z_coordinates, const, 'pos');
A_estimates = [];
B_estimates = [];
C_estimates = [];
D_estimates = [];
for i = 1 : size(all_tf_j_pos, 2)
    [A_est, B_est, C_est, D_est] = dra(all_tf_j_pos(:, i), res0, D, sampling_freq, T_len, const);
    A_estimates = [A_estimates; A_est];
    B_estimates = [B_estimates; B_est];
    C_estimates = [C_estimates; C_est];
    D_estimates = [D_estimates, D_est];
end


steps = 10000;
li_flux_pos = zeros(steps, size(all_tf_j_pos, 2));
X = zeros(size(A_estimates, 2), size(all_tf_j_pos, 2));
for step = 1 : steps
    for i = 1 : size(all_tf_j_pos, 2)
        A = A_estimates(i * size(A_estimates, 2) - size(A_estimates, 2) + 1: i * size(A_estimates, 2), 1 : size(A_estimates, 2));
        B = B_estimates(i * size(A_estimates, 2) - size(A_estimates, 2) + 1: i * size(A_estimates, 2), 1);
        C = C_estimates(i, 1 : size(A_estimates, 2));
        D = D_estimates(1, i);
        X(:, i) = A * X(:, i) + B;
        Y = C * X(:, i) + D;
        li_flux_pos(step, i) = Y;
    end
end

t = 1 : steps;
subplot(2, 1, 1);
for i = 1 : size(li_flux_pos, 2)
    plot(t, li_flux_pos(:, i))
    text(max(t), max(li_flux_pos(:, i)), "Pos" + num2str(z_coordinates(i)));
    hold on;
end
grid on;

[all_tf_j_neg, res0, D, sampling_freq, T_len] = tf_j(cse_neg, z_coordinates, const, 'neg');
A_estimates = [];
B_estimates = [];
C_estimates = [];
D_estimates = [];
for i = 1 : size(all_tf_j_neg, 2)
    [A_est, B_est, C_est, D_est] = dra(all_tf_j_neg(:, i), res0, D, sampling_freq, T_len, const);
    A_estimates = [A_estimates; A_est];
    B_estimates = [B_estimates; B_est];
    C_estimates = [C_estimates; C_est];
    D_estimates = [D_estimates, D_est];
end


steps = 10000;
li_flux_neg = zeros(steps, size(all_tf_j_neg, 2));
X = zeros(size(A_estimates, 2), size(all_tf_j_neg, 2));
for step = 1 : steps
    for i = 1 : size(all_tf_j_neg, 2)
        A = A_estimates(i * size(A_estimates, 2) - size(A_estimates, 2) + 1: i * size(A_estimates, 2), 1 : size(A_estimates, 2));
        B = B_estimates(i * size(A_estimates, 2) - size(A_estimates, 2) + 1: i * size(A_estimates, 2), 1);
        C = C_estimates(i, 1 : size(A_estimates, 2));
        D = D_estimates(1, i);
        X(:, i) = A * X(:, i) + B;
        Y = C * X(:, i) + D;
        li_flux_neg(step, i) = Y;
    end
end

t = 1 : steps;
subplot(2, 1, 2);
for i = 1 : size(li_flux_neg, 2)
    plot(t, li_flux_neg(:, i))
    text(max(t), max(li_flux_neg(:, i)), "Neg" + num2str(z_coordinates(i)));
    hold on;
end
grid on;






% steps = 10000;
% li_flux_vector_pos = zeros(steps, 1);
% X_li_flux = zeros(size(A_est, 2), 1);
% for current_step = 1 : steps
%     X_li_flux = A_est * X_li_flux + B_est;
%     Y_li_flux = C_est * X_li_flux + D_est;
%     if current_step == 1
%         li_flux_vector_pos(current_step, 1) = 0;
%     else
%         li_flux_vector_pos(current_step, 1) = li_flux_vector_pos(current_step - 1, 1) + Y_li_flux;
%     end
% end
% t_vector = 1 : size(li_flux_vector_pos, 1);
% plot(t_vector, li_flux_vector_pos)
% hold on;
% text(max(t_vector), min(li_flux_vector_pos), 'pos')
% 
% [tf_j_neg, res0, D, sampling_freq, T_len] = tf_j(cse_neg, z_coordinates, const, 'neg');
% [A_est, B_est, C_est, D_est] = dra(tf_j_neg, res0, D, sampling_freq, T_len, const);
% 
% 
% li_flux_vector_neg = zeros(steps, 1);
% X_li_flux = zeros(size(A_est, 2), 1);
% for current_step = 1 : steps
%     X_li_flux = A_est * X_li_flux + B_est;
%     Y_li_flux = C_est * X_li_flux + D_est;
%     if current_step == 1
%         li_flux_vector_neg(current_step, 1) = 0;
%     else
%         li_flux_vector_neg(current_step, 1) = li_flux_vector_neg(current_step - 1, 1) + Y_li_flux;
%     end
% end
% t_vector = 1 : size(li_flux_vector_neg, 1);
% plot(t_vector, li_flux_vector_neg)
% text(max(t_vector), max(li_flux_vector_neg), 'neg')
% grid on;




fprintf("\n\n");




% Initial state variables.
% cse_neg = 0;
% cs0_neg = 0;

% Estimate ss for different tf's and states.
% [j_neg_ss] = estimate_li_flux_neg_ss(cse_neg, const);
% 
% [A_cse_neg, B_cse_neg, C_cse_neg, D_cse_neg] = estimate_cse_neg_ss(const);
% sys_cse_neg = ss(A_cse_neg, B_cse_neg, C_cse_neg, D_cse_neg);
% 
% li_flux_neg_vector = [0.01; zeros(size(load_cycle, 1) - 1, 1)];
% 
% X_cse_neg = zeros(size(A_cse_neg, 1), 1);
% X_li_flux_neg = zeros(size(A_li_neg, 1), 1);
% for current_step = 1 : size(load_cycle, 1)
%     current_load = load_cycle(current_step);
%     X_li_flux_neg = A_li_neg * X_li_flux_neg + B_li_neg * current_load;
%     Y_li_flux_neg = C_li_neg * X_li_flux_neg + D_li_neg * current_load;
%     if current_step == 1
%         flux_neg_vector(current_step, 1) = 0;
%     else
%         li_flux_neg_vector(current_step, 1) = li_flux_neg_vector(current_step - 1, 1) + Y_li_flux_neg;
%     end
% 
%     % cse_neg state space.
%     X_cse_neg = A_cse_neg * X_cse_neg + B_cse_neg * current_load;
%     delta_cse_neg = C_cse_neg * X_cse_neg + D_cse_neg * current_load;
%     cse_neg = cse_neg + delta_cse_neg / size(load_cycle, 1);
% 
%     % Calculate the ocv_neg derivative.
%     ocv_derivative_neg = calculate_ocv_derivative_neg(cse_neg, const);
% end
% t_vector = 1 : size(li_flux_neg_vector, 1);
% % plot(t_vector, li_flux_neg_vector)
% 
% fprintf("\n\n");
