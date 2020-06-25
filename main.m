clearvars;
run('parameters.m');

cse_neg = 2000;
z_coordinates = [0.0 .3 0.6 1.0];
[tf_j_neg, res0, D, sampling_freq, T_len] = tf_j(cse_neg, z_coordinates, const, 'neg');
[tf_j_pos, res0, D, sampling_freq, T_len] = tf_j(cse_neg, z_coordinates, const, 'pos');
% t_vector = 1 : size(tf_j, 2)
% plot(t_vector, tf_j_neg)

[A_est, B_est, C_est, D_est] = dra(tf_j_pos, res0, sampling_freq, T_len, const);
% sys = ss(A_est, B_est, C_est, D_est, 0.01);
% impulse(sys)


steps = 1000;
li_flux_vector = zeros(steps, 1);
X_li_flux = zeros(size(A_est, 2), 1);
for current_step = 1 : steps
    X_li_flux = A_est * X_li_flux + B_est;
    Y_li_flux = C_est * X_li_flux + D_est;
    if current_step == 1
        li_flux_vector(current_step, 1) = 0;
    else
        li_flux_vector(current_step, 1) = li_flux_vector(current_step - 1, 1) + Y_li_flux;
    end
end
t_vector = 1 : size(li_flux_vector, 1);
plot(t_vector, li_flux_vector)
hold on;

[A_est, B_est, C_est, D_est] = dra(tf_j_neg, res0, sampling_freq, T_len, const);

li_flux_vector = zeros(steps, 1);
X_li_flux = zeros(size(A_est, 2), 1);
for current_step = 1 : steps
    X_li_flux = A_est * X_li_flux + B_est;
    Y_li_flux = C_est * X_li_flux + D_est;
    if current_step == 1
        li_flux_vector(current_step, 1) = 0;
    else
        li_flux_vector(current_step, 1) = li_flux_vector(current_step - 1, 1) + Y_li_flux;
    end
end
t_vector = 1 : size(li_flux_vector, 1);
plot(t_vector, li_flux_vector)


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
