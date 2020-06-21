clearvars;
run('parameters.m');

% Initial state variables.
cse_neg = 1000;
cs0_neg = 1000;

% Load cycle.
load_cycle_len = 1000;
i_app = 1;
load_cycle = i_app * ones(load_cycle_len, 1);

% Estimate ss for different tf's and states.
[A_li_neg, B_li_neg, C_li_neg, D_li_neg] = estimate_li_flux_neg_ss(cse_neg, solid_max_concentration_neg, sigma_eff_neg, k_eff_neg, L_neg, asymmetric_charge_transfer_neg, resistivity_neg, ionic_conductivity_neg, R_solid_electrolyte_neg, R_neg, diffusivity_neg, F, A_neg);

disp(A_li_neg)
% [A_cse_neg, B_cse_neg, C_cse_neg, D_cse_neg] = estimate_cse_neg_ss(diffusivity_neg, radius_neg);

% X = zeros(size(A_cse_neg, 1), 1);
% for current_step = 1 : size(load_cycle, 1)
%     current_load = load_cycle(current_step);
%     % cse_neg state space.
%     X = A_cse_neg * X + B_cse_neg * current_load;
%     delta_cse_neg = C_cse_neg * X + D_cse_neg * load_cycle(current_step);
%     cse_neg = cse_neg + delta_cse_neg / size(load_cycle, 1);
% 
%     % Calculate the ocv_neg derivative.
%     ocv_derivative_neg = calculate_ocv_derivative_neg(cse_neg, solid_max_concentration_neg);
% end
% fprintf("\n\n\n");
