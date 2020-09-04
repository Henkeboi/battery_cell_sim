clearvars;
disp("Loading data.")
run('parameters.m');
run('read_load_cycle.m')

%%% Calculations for the negative eletrode at z = 0.
cs0_neg = const.cs0_neg;
cs0_pos = const.cs0_pos;
cs_max_neg = const.solid_max_c_neg;
cs_max_pos = const.solid_max_c_pos;
ce0_neg = const.ce0_neg;
ce0_pos = const.ce0_pos;


%f = 200;
%T_len = 20;
%Ts = 0.5;
%z = 0.0;
%
%[tf_delta_cse, res0, D] = tf_delta_cse(90000, z, T_len, f, 'pos', const);
%[A, B, C, D] = dra(tf_delta_cse, res0, f, T_len, Ts, const);
%S = ss(A, B, C, D, Ts);
%
%disp("Is stable: " + isstable(S))
%%step(S)
%
%P = 2;
%S = 144;
%time_acc = 0;
%time = zeros(size(load_cycle, 1), 1); 
%measured_voltage = zeros(size(load_cycle, 1), 1);
%X = zeros(5, 1);
%Y = zeros(size(load_cycle, 1), 1);
%cse0 = 90000;
%a = zeros(size(load_cycle, 1), 1);
%for i = 1 : size(load_cycle, 1)
%    % Find U.
%    delta_time = load_cycle(i, 1);
%    current = load_cycle(i, 2);
%    measured_voltage(i) = load_cycle(i, 3) / S;
%    U = current * delta_time / P;
%    X = A * X + B * U;
%    Y(i) = C * X + D * U;  
%
%    time(i) = time_acc;
%    time_acc = time(i) + delta_time;
%    cse0 = cse0 - Y(i);
%    a(i) = cse0;
%end
%plot(time, a)
%return

disp("Running blended Lithium flux dra.")
j_neg_sampling_f = 400;
j_neg_T_len = 10;
j_neg_Ts = 0.1;
j_blender_neg = Blender(0.1, @tf_j, j_neg_Ts, 0, 'neg', const);
j_blender_neg.create_models(j_neg_T_len, j_neg_sampling_f);
j_blender_neg.sort();

j_pos_sampling_f = 400;
j_pos_T_len = 10;
j_pos_Ts = 0.1;
j_blender_pos = Blender(0.1, @tf_j, j_pos_Ts, 0, 'pos', const);
j_blender_pos.create_models(j_pos_T_len, j_pos_sampling_f);
j_blender_pos.sort();

disp("Running blended cse dra.")
cse_neg_sampling_f = 400;
cse_neg_T_len = 10;
cse_neg_Ts = 0.1;
cse_blender_neg = Blender(0.1, @tf_cse, cse_neg_Ts, 0.0, 'neg', const);
cse_blender_neg.create_models(cse_neg_T_len, cse_neg_sampling_f);
cse_blender_neg.sort();

cse_pos_sampling_f = 400;
cse_pos_T_len = 10;
cse_pos_Ts = 0.1;
cse_blender_pos = Blender(0.1, @tf_cse, cse_pos_Ts, 0.1, 'pos', const);
cse_blender_pos.create_models(cse_pos_T_len, cse_pos_sampling_f);
cse_blender_pos.sort();

%T_len_cse_neg = 20;
%Ts_cse_neg = 0.1;
%z_cse_neg = 0;
%[tf_delta_cse_neg, res0_cse_neg, D_cse_neg] = tf_delta_cse(10000, z_cse_neg, T_len_cse_neg, f_cse_neg, 'neg', const);
%[A_cse_neg, B_cse_neg, C_cse_neg, D_cse_neg] = dra(tf_delta_cse_neg, res0_cse_neg, f_cse_neg, T_len_cse_neg, Ts_cse_neg, const);
%
%f_cse_pos = 200;
%T_len_cse_pos = 10;
%Ts_cse_pos = 10;
%z_cse_pos = 0;
%[tf_delta_cse_pos, res0_cse_pos, D_cse_pos] = tf_delta_cse(10000, z_cse_pos, T_len_cse_pos, f_cse_pos, 'pos', const);
%[A_cse_pos, B_cse_pos, C_cse_pos, D_cse_pos] = dra(tf_delta_cse_pos, res0_cse_pos, f_cse_pos, T_len_cse_pos, Ts_cse_pos, const);

disp("Running blended ce dra.")
ce_neg_sampling = 200;
ce_neg_T_len = 10;
ce_neg_Ts = 0.1;
x = 0; 
ce_blender_neg = Blender(0.1, @tf_ce, ce_neg_Ts, x, 'neg', const);
ce_blender_neg.create_models(ce_neg_T_len, ce_neg_sampling);
ce_blender_neg.sort();

ce_pos_sampling = 200;
ce_pos_T_len = 10;
ce_pos_Ts = 0.1;
x = const.L_neg + const.L_sep + const.L_pos; 
ce_blender_pos = Blender(0.1, @tf_ce, ce_pos_Ts, x, 'pos', const);
ce_blender_pos.create_models(ce_pos_T_len, ce_pos_sampling);
ce_blender_pos.sort();

disp("Running blended pote dra.")
pote_L_tot_sampling = 200;
pote_L_tot_T_len = 10;
pote_L_tot_Ts = 0.1;
x = const.L_neg + const.L_sep + const.L_pos; 
pote_blender_L_tot = Blender(0.1, @tf_pote, pote_L_tot_Ts, x, 'neg', const);
pote_blender_L_tot.create_models(pote_L_tot_T_len, pote_L_tot_sampling);
pote_blender_L_tot.sort();

disp("Simulating.")
z_neg = zeros(size(load_cycle, 1), 1);
z_pos = zeros(size(load_cycle, 1), 1);
j_neg = zeros(size(load_cycle, 1), 1);
j_pos = zeros(size(load_cycle, 1), 1);
pots_neg = zeros(size(load_cycle, 1), 1);
pots_pos = zeros(size(load_cycle, 1), 1);
cse_neg = zeros(size(load_cycle, 1), 1);
cse_pos = zeros(size(load_cycle, 1), 1);
ce_neg = zeros(size(load_cycle, 1), 1);
ce_pos = zeros(size(load_cycle, 1), 1);
potse_neg = zeros(size(load_cycle, 1), 1);
potse_pos = zeros(size(load_cycle, 1), 1);
pote_L_tot_1 = zeros(size(load_cycle, 1), 1);
pote_L_tot_2 = zeros(size(load_cycle, 1), 1);
v = zeros(size(load_cycle, 1), 1);
measured_voltage = zeros(size(load_cycle, 1), 1);
time = zeros(size(load_cycle, 1), 1);
time_acc = 0;
S = 144;
P = 2;

X_cse_neg = zeros(5, 1);
X_cse_pos = zeros(5, 1);
SOC_pos = 1;
SOC_neg = 1;

for i = 1 : size(load_cycle, 1)
    % Find U.
    delta_time = load_cycle(i, 1);
    current = load_cycle(i, 2);
    measured_voltage(i) = load_cycle(i, 3) / S;
    U = current * delta_time / P;

    [j_neg_X, j_neg_Y, j_neg_integrator_index] = j_blender_neg.step(U, SOC_neg);
    j_neg(i) = j_neg_Y;
    [j_pos_X, j_pos_Y, j_pos_integrator_index] = j_blender_pos.step(U, SOC_pos);
    j_pos(i) = j_pos_Y;

    [cse_neg_X, cse_neg_Y, cse_neg_integrator_index] = cse_blender_neg.step(U, SOC_neg);
    %cse_neg(i) = cse_neg_Y + cs0_neg; % Debias
    X_cse_neg_test = A_cse_neg * X_cse_neg + B_cse_neg * U; % Test
    cse_neg_test(i) = C_cse_neg * X_cse_neg + D_cse_neg * U; % Test
    cs_avg_neg(i) = cs0_neg - cse_neg_X(cse_neg_integrator_index);
    z_neg(i) = (cs_avg_neg(i) - const.solid_max_c_neg * const.x0_neg) / (const.solid_max_c_neg * const.x100_neg - const.solid_max_c_neg * const.x0_neg);
    SOC_neg = z_neg(i);

    [cse_pos_X, cse_pos_Y, cse_pos_integrator_index] = cse_blender_pos.step(U, SOC_pos);
    %cse_pos(i) = cse_pos_Y + cs0_pos; % Debias
    X_cse_pos_test = A_cse_pos * X_cse_pos + B_cse_pos; % Test
    cse_pos_test(i) = C_cse_pos * X_cse_pos + D_cse_pos; % Test
    cs_avg_pos(i) = cs0_pos + cse_pos_X(cse_pos_integrator_index);
    z_pos(i) = 1 - (cs_avg_pos(i) - const.solid_max_c_pos * const.x0_pos) / (const.solid_max_c_pos * const.x100_pos - const.solid_max_c_pos * const.x0_pos);
    SOC_pos = z_pos(i);

    [ce_neg_X, ce_neg_Y, ce_neg_integrator_index] = ce_blender_neg.step(U, SOC_neg);
    ce_neg(i) = ce_neg_Y + ce0_neg;
    [ce_pos_X, ce_pos_Y, ce_pos_integrator_index] = ce_blender_pos.step(U, SOC_pos);
    ce_pos(i) = ce_pos_Y + ce0_pos;

    [pote_L_tot_X, pote_L_tot_Y, pote_L_tot_integrator_index] = pote_blender_L_tot.step(U, SOC_pos);
    pote_L_tot_1(i) = pote_L_tot_Y;
    pote_L_tot_2(i) = 2 * const.R * const.temp * (1 - const.transfer_number_pos) * (log(ce_pos(i)) - log(ce0_pos)) / const.F;

    %v(i) = calculate_voltage(j_neg(i), j_pos(i), cse_neg(i), cse_pos(i), ce_neg(i), ce_pos(i), pote_L_tot_1(i), pote_L_tot_2(i), const);
    v(i) = calculate_voltage(j_neg(i), j_pos(i), cse_neg_test(i), cse_pos_test(i), ce_neg(i), ce_pos(i), pote_L_tot_1(i), pote_L_tot_2(i), const); % Test

    time(i) = time_acc;
    time_acc = time(i) + delta_time;
end

f1 = figure;
plot(time, v);
hold on;
plot(time, measured_voltage);
title("Voltage")
xlabel("Time")
ylabel("[V]")
grid on;

%f2 = figure;
%plot(time, z_neg);
%hold on;
%plot(time, z_pos, 'r');
%title("SOC")
%xlabel("Time")
%ylabel("SOC")
%grid on;

%f2 = figure;
%plot(time, cse_neg);
%hold on;
%plot(time, cse_pos, 'r');
%title("Average surface concentration")
%xlabel("Time")
%ylabel("Lithium concentration [mol / m^3]")
%grid on;

%f3 = figure;
%plot(time, potse_neg);
%hold on;
%plot(time, potse_pos, 'r');
%title("potse")
%xlabel("Time")
%ylabel("[V]")
%grid on;

%f5 = figure;
%plot(time, ce_neg);
%hold on;
%plot(time, ce_pos, 'r');
%title("ce")
%xlabel("Time")
%ylabel("Lithium concentration [mol / m^3]")
%grid on;

%f6 = figure;
%plot(time, pote_L_tot_1);
%hold on;
%plot(time, pote_L_tot_2, 'r');
%title("Pote")
%xlabel("Time")
%ylabel("[V]")
%grid on;

%f6 = figure;
%plot(time, j_neg);
%hold on;
%plot(time, j_pos, 'r');
%title("j")
%xlabel("Time")
%ylabel("[V]")
%grid on;
