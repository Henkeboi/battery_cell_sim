clearvars;
run('parameters.m');
[A_cse_neg, B_cse_neg, C_cse_neg, D_cse_neg] = estimate_cse_neg_ss(diffusivity_neg, radius_neg);

load_cycle_len = 10;
i_app = 1;
load_cycle = i_app * ones(load_cycle_len, 1);
X = zeros(size(A_cse_neg, 1), 1);
for current_step = 1 : size(load_cycle, 1)
    current_load = load_cycle(current_step);
    % cse_neg state space.
    X = A_cse_neg * X / current_load + B_cse_neg * current_load;
    delta_cse_neg = C_cse_neg * X + D_cse_neg * load_cycle(current_step);
    cse_neg = cse_neg + delta_cse_neg / size(load_cycle, 1);

    % Calculate the ocv_neg derivative.
    ocv_derivative_neg = calculate_ocv_derivative_neg(cse_neg, solid_max_concentration_neg);
    disp(ocv_derivative_neg);
    disp(cse_neg);
end
fprintf("\n\n\n");
