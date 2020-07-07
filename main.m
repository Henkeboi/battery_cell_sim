clearvars;
disp("Loading data.")
run('parameters.m');
run('read_load_cycle.m')

%% Calculations for the negative eletrode at z = 0.
cs0_neg = const.cse0_neg;
cs0_pos = const.cse0_pos;
disp("Running blended Lithium surface concentration dra.")

cse_neg_sampling_f = 400;
cse_neg_T_len = 100;
cse_neg_blender = Blender(0.01, @tf_cse, [0]);
[all_cse_neg_A, all_cse_neg_B, all_cse_neg_C, all_cse_neg_D, cse_neg_integrator_index, cse_neg_Ts] = cse_neg_blender.create_models(cse_neg_T_len, cse_neg_sampling_f, 'neg', const);
[cse_neg_A, cse_neg_B, cse_neg_C, cse_neg_D] = cse_neg_blender.blend_model(all_cse_neg_A, all_cse_neg_B, all_cse_neg_C, all_cse_neg_D, calculate_SOC(cs0_neg, 0, 'neg', const));

cse_pos_sampling_f = 400;
cse_pos_T_len = 5;
cse_pos_blender = Blender(0.01, @tf_cse, [0]);
[all_cse_pos_A, all_cse_pos_B, all_cse_pos_C, all_cse_pos_D, cse_pos_integrator_index, cse_pos_Ts] = cse_pos_blender.create_models(cse_pos_T_len, cse_pos_sampling_f, 'pos', const);
[cse_pos_A, cse_pos_B, cse_pos_C, cse_pos_D] = cse_pos_blender.blend_model(all_cse_pos_A, all_cse_pos_B, all_cse_pos_C, all_cse_pos_D, calculate_SOC(cs0_pos, 0, 'pos', const));

% disp("Running blended Lithium flux dra.")
% j_blender = Blender(0.1, @tf_cse, [0]);
% [all_j_A, all_j_B, all_j_C, all_j_D, j_integrator_index, Ts] = j_blender.create_models('neg', const);
% [j_A, j_B, j_C, j_D] = j_blender.blend_model(all_j_A, all_j_B, all_j_C, all_j_D, calculate_SOC(cs0, 0, 'neg', const));
% 
% disp("Running blended pots dra.")
% pots_blender = Blender(0.1, @tf_pots, [0]);
% [all_pots_A, all_pots_B, all_pots_C, all_pots_D, pots_integrator_index, Ts] = pots_blender.create_models('neg', const);
% [pots_A, pots_B, pots_C, pots_D] = pots_blender.blend_model(all_pots_A, all_pots_B, all_pots_C, all_pots_D, calculate_SOC(cs0, 0, 'neg', const));
% 
% disp("Running blended potse dra.")
% potse_blender = Blender(0.1, @tf_potse, [0]);
% [all_potse_A, all_potse_B, all_potse_C, all_potse_D, potse_integrator_index, Ts] = potse_blender.create_models('neg', const);
% [potse_A, potse_B, potse_C, potse_D] = pots_blender.blend_model(all_potse_A, all_potse_B, all_potse_C, all_potse_D, calculate_SOC(cs0, 0, 'neg', const));
% 
% disp("Running unblended ce dra.")
% locs = [2e-4];
% [tf_ce, res0, D, sampling_freq, T_len] = tf_ce(0.6, 0.6, 5000, 5000, @calculate_ocv_derivative_neg, @calculate_ocv_derivative_pos, locs, const);
% [ce_A, ce_B, ce_C, ce_D, ce_Ts] = dra(tf_ce, res0, D, sampling_freq, T_len, const);
% [ce_A, ce_B, ce_C, ce_D, ce_integrator_index] = multi_dra(ce_A, ce_B, ce_C, ce_D, ce_Ts, res0);

disp("Simulating.")
cse_neg_X = zeros(size(cse_neg_A, 1), 1);
cse_pos_X = zeros(size(cse_pos_A, 1), 1);
% j_X = zeros(size(j_A, 1), 1);
% pots_X = zeros(size(pots_A, 1), 1);
% potse_X = zeros(size(potse_A, 1), 1);
% ce_X = zeros(size(ce_A, 1), 1);

z_neg = zeros(size(load_cycle, 1), 1);
z_pos = zeros(size(load_cycle, 1), 1);
cse_neg = zeros(size(load_cycle, 1), 1);
cse_pos = zeros(size(load_cycle, 1), 1);
% j = zeros(size(load_cycle, 1), 1);
% pots = zeros(size(load_cycle, 1), 1);
% potse = zeros(size(load_cycle, 1), 1);
% ce = zeros(size(load_cycle, 1), 1);
time = zeros(size(load_cycle, 1), 1);
time_acc = 0;

% TODO: Debias the variables.
for i = 1 : size(load_cycle, 1)
    % Find U.
    delta_time = load_cycle(i, 1);
    current = load_cycle(i, 2);
    P = 2;
    U = current * delta_time / cse_neg_Ts / P;

    % Find SOC and lithium surface concentration.
    cse_neg_X = cse_neg_A * cse_neg_X + cse_neg_B * U;
    cse_neg(i) = cs0_neg + cse_neg_X(cse_neg_integrator_index);
    z_neg(i) = calculate_SOC(cs0_neg, cse_neg_X(cse_neg_integrator_index), 'neg', const); 

    cse_pos_X = cse_pos_A * cse_pos_X + cse_pos_B * U;
    cse_pos(i) = cs0_pos + cse_pos_X(cse_pos_integrator_index);
    z_pos(i) = calculate_SOC(cs0_pos, cse_pos_X(cse_pos_integrator_index), 'pos', const); 

    % Find li flux
    % j_X = j_A * j_X + j_B * U;
    % j(i) = j_X(j_integrator_index);

    % Find pots
    % pots_X = pots_A * pots_X + pots_C * U;
    % pots(i) = pots_X(pots_integrator_index);

    % Find potse
    % potse_X = potse_A * potse_X + potse_C * U;
    % potse(i) = potse_X(potse_integrator_index);

    % Find ce
    % ce_X = ce_A * ce_X + ce_B * U;
    % ce_Y = ce_C * ce_X; 
    % ce(i) = ce_Y + const.ce0;

    % Find voltage
    % disp(cs0 - cse_X(cse_integrator_index) / const.porosity_solid_neg / const.A_neg / const.L_neg / const.F)
    % param.cse_neg = Y_cse;
    % v = calculate_voltage(param, const);

    % Find next step state space.
    [cse_neg_A, cse_neg_B, cse_neg_C, cse_neg_D] = cse_neg_blender.blend_model(all_cse_neg_A, all_cse_neg_B, all_cse_neg_C, all_cse_neg_D, z_neg(i));
    [cse_pos_A, cse_pos_B, cse_pos_C, cse_pos_D] = cse_pos_blender.blend_model(all_cse_pos_A, all_cse_pos_B, all_cse_pos_C, all_cse_pos_D, z_pos(i));

    % [j_A, j_B, j_C, j_D] = j_blender.blend_model(all_j_A, all_j_B, all_j_C, all_j_D, z(i));
    % [pots_A, pots_B, pots_C, pots_D] = pots_blender.blend_model(all_pots_A, all_pots_B, all_pots_C, all_pots_D, z(i));
    % [potse_A, potse_B, potse_C, potse_D] = potse_blender.blend_model(all_potse_A, all_potse_B, all_potse_C, all_potse_D, z(i));
    time(i) = time_acc;
    time_acc = time(i) + delta_time;
end
disp("Plotting results.")

f1 = figure;
plot(time, z_neg);
hold on;
plot(time, z_pos, 'r');
title("SOC")
xlabel("Time")
ylabel("SOC")
grid on;

f2 = figure;
plot(time, cse_neg);
hold on;
plot(time, cse_pos, 'r');
title("Average concentration")
xlabel("Time")
ylabel("Lithium concentration [mol / m^3]")
grid on;


% 
% f4 = figure;
% plot(time, cse_pos);
% title("Average concentration at the positive electrode")
% xlabel("Time")
% ylabel("Lithium concentration [mol / m^3]")
% grid on;

% f3 = figure;
% plot(time, z);
% title("SOC")
% xlabel("Time")
% ylabel("SOC")

% f4 = figure;
% plot(time, j);
% title("Li flux at z = 0")
% xlabel("Time")
% ylabel("Li flux")
% 
% f5 = figure;
% plot(time, pots);
% title("Solid potential at z = 0")
% xlabel("Time")
% ylabel("Solid potential")
% 
% f6 = figure;
% plot(time, potse);
% title("Solid electrolyte potential at z = 0")
% xlabel("Time")
% ylabel("Solid potential")
% 
% f7 = figure;
% plot(time, ce);
% title("Lithium concentration in the electrolyte at x = 0")
% xlabel("Time")
% ylabel("Lithium concentration")
