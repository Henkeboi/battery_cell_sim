clearvars;
disp("Loading data.")
run('parameters.m');
run('read_load_cycle.m')

%% Calculations for the negative eletrode at z = 0.
cs0_neg = const.solid_max_c_neg * const.x100_neg;
cs0_pos = const.solid_max_c_pos * const.x0_pos;
cs_max_neg = const.solid_max_c_neg;
cs_max_pos = const.solid_max_c_pos;

% cse = 20000;
% [tf_potse, res0]  = tf_pots(cse, [0], 100, 200, 'neg', const);
% [A, B, C, D, T_shifted] = dra(tf_potse, res0, D, 200, 100, const);
% S = ss(A, B, C, D, T_shifted);

disp("Running blended Lithium surface concentration dra.")
cse_neg_sampling_f = 200;
cse_neg_T_len = 10;
cse_neg_Ts = 0.5;
cse_neg_blender = Blender(0.1, @tf_cse, cse_neg_Ts, [0], 'neg', const);
cse_neg_blender.create_models(cse_neg_T_len, cse_neg_sampling_f);
cse_neg_blender.sort();

cse_pos_sampling_f = 200;
cse_pos_T_len = 10;
cse_pos_Ts = 0.5;
cse_pos_blender = Blender(0.1, @tf_cse, cse_pos_Ts, [0], 'pos', const);
cse_pos_blender.create_models(cse_pos_T_len, cse_pos_sampling_f);
cse_pos_blender.sort();
 
disp("Running blended lithium flux dra.")
j_neg_sampling = 200;
j_neg_T_len = 10;
j_neg_Ts = 0.1;
j_neg_blender = Blender(0.1, @tf_j, j_neg_Ts, [0], 'neg', const);
j_neg_blender.create_models(j_neg_T_len, j_neg_sampling);
j_neg_blender.sort();

j_pos_sampling = 200;
j_pos_T_len = 10;
j_pos_Ts = 0.1;
j_pos_blender = Blender(0.1, @tf_j, j_pos_Ts, [0], 'pos', const);
j_pos_blender.create_models(j_pos_T_len, j_pos_sampling);
j_pos_blender.sort();

disp("Running blended potse dra.")
potse_neg_sampling = 400;
potse_neg_T_len = 100;
potse_neg_Ts = 0.1;
potse_neg_blender = Blender(0.1, @tf_potse, potse_neg_Ts, [0], 'neg', const);
potse_neg_blender.create_models(potse_neg_T_len, potse_neg_sampling);
potse_neg_blender.sort();

potse_pos_sampling = 200;
potse_pos_T_len = 10;
potse_pos_Ts = 0.5;
potse_pos_blender = Blender(0.1, @tf_potse, potse_pos_Ts, [0], 'pos', const);
potse_pos_blender.create_models(potse_pos_T_len, potse_pos_sampling);
potse_pos_blender.sort();

disp("Running blended pots dra.")
pots_neg_sampling = 200;
pots_neg_T_len = 10;
potse_neg_Ts = 0.5;
pots_neg_blender = Blender(0.1, @tf_pots, potse_neg_Ts, [0], 'neg', const);
pots_neg_blender.create_models(pots_neg_T_len, pots_neg_sampling);
pots_neg_blender.sort();

pots_pos_sampling = 200;
pots_pos_T_len = 10;
potse_pos_Ts = 0.5;
pots_pos_blender = Blender(0.1, @tf_pots, potse_pos_Ts, [0], 'pos', const);
pots_pos_blender.create_models(pots_pos_T_len, pots_pos_sampling);
pots_pos_blender.sort();

disp("Simulating.")
z_neg = zeros(size(load_cycle, 1), 1);
z_pos = zeros(size(load_cycle, 1), 1);
cse_neg = zeros(size(load_cycle, 1), 1);
cse_pos = zeros(size(load_cycle, 1), 1);
j_neg = zeros(size(load_cycle, 1), 1);
j_pos = zeros(size(load_cycle, 1), 1);
potse_neg = zeros(size(load_cycle, 1), 1);
potse_pos = zeros(size(load_cycle, 1), 1);
pots_neg = zeros(size(load_cycle, 1), 1);
pots_pos = zeros(size(load_cycle, 1), 1);
time = zeros(size(load_cycle, 1), 1);
v = zeros(size(load_cycle, 1), 1);
time_acc = 0;

SOC_pos = 1;
SOC_neg = 1;
for i = 1 : size(load_cycle, 1)
    % Find U.
    delta_time = load_cycle(i, 1);
    current = load_cycle(i, 2);
    P = 2;
    U = current * delta_time / P;

    [cse_neg_X, cse_neg_Y, cse_neg_integrator_index] = cse_neg_blender.step(U, SOC_neg);
    cse_neg(i) = cse_neg_Y;
    z_neg(i) = calculate_SOC(cs0_neg, cse_neg_X(cse_neg_integrator_index), 'neg', const); 
    SOC_neg = z_neg(i);
    cs_avg_neg = cs0_neg - cse_neg_X(cse_neg_integrator_index);

    [cse_pos_X, cse_pos_Y, cse_pos_integrator_index] = cse_pos_blender.step(U, SOC_pos);
    cse_pos(i) = cse_pos_Y;
    z_pos(i) = calculate_SOC(cs0_pos, cse_pos_X(cse_pos_integrator_index), 'pos', const); 
    SOC_pos = z_pos(i);
    cs_avg_pos = cs0_pos + cse_pos_X(cse_pos_integrator_index);

    [j_neg_X, j_neg_Y, j_neg_integrator_index] = j_neg_blender.step(U, SOC_neg);
    j_neg(i) = j_neg_Y;
    [j_pos_X, j_pos_Y, j_pos_integrator_index] = j_pos_blender.step(U, SOC_pos);
    j_pos(i) = j_pos_Y;
 
    %[potse_neg_X, potse_neg_Y, potse_neg_integrator_index] = potse_neg_blender.step(U, SOC_neg);
    %potse_neg(i) = potse_neg_Y;
    %[potse_pos_X, potse_pos_Y, potse_pos_integrator_index] = potse_pos_blender.step(U, SOC_pos);
    %potse_pos(i) = potse_pos_Y + potse_pos_X(potse_pos_integrator_index); 

    %[pots_neg_X, pots_neg_Y, pots_neg_integrator_index] = pots_neg_blender.step(U, SOC_neg);
    %pots_neg(i) = pots_neg_Y;
    if i == 10
        % [A, B, C, D, Ts] = pots_neg_blender.blend_model(SOC_neg);
        % S = ss(A, B, C, D, potse_neg_Ts);
        % [A, B, C, D, Ts] = potse_pos_blender.blend_model(SOC_pos);
        % S = ss(A, B, C, D, -1);
        % pzmap(S);
        % bode(S)
        % return;
        % [y, k] = impulse(S, 0:1000);
        % stem(y, k, 'filled')
        %return;
        %pzmap(S)
        % bode(S)
        %return
    end
    % %[pots_pos_X, pots_pos_Y, pots_pos_integrator_index] = pots_pos_blender.step(U, SOC_pos);
    % pots_pos(i) = pots_pos_Y + pots_pos_X(pots_pos_integrator_index); 

    % pote1 = 0;
    % pote2 = 0;
    % ce_neg = cse_neg(i);
    % ce_pos = cse_pos(i);
    % v(i) = calculate_voltage(cse_neg(i), cse_pos(i), j_neg(i), j_pos(i), pote1, pote2, ce_neg, ce_pos, const);

    time(i) = time_acc;
    time_acc = time(i) + delta_time;

    if isnan(z_neg(i))
        disp("Nan z_neg")
        return;
    end
    if isnan(z_pos(i))
        disp("Nan z_pos")
        return;
    end
    if isnan(cse_neg(i))
        disp("Nan cse_neg")
        break;
    end
    if isnan(cse_pos(i))
        disp(i);
        disp("NAN cse_pos")
        break;
    end
    if isnan(j_neg(i))
        disp("NAN j_neg")
        break;
    end
    if isnan(j_pos(i))
        disp("NAN j_neg")
        break;
    end
    if isnan(potse_neg(i))
        disp("NAN potse_neg")
        return;
    end
    if isnan(potse_pos(i))
        disp("NAN potse_pos")
        return;
    end
    if isnan(pots_neg(i))
        disp("NAN pots_neg")
        return;
    end
    if isnan(pots_pos(i))
        disp("NAN pots_pos")
        return;
    end
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

% f2 = figure;
% plot(time, cse_neg);
% hold on;
% plot(time, cse_pos, 'r');
% title("Average surface concentration")
% xlabel("Time")
% ylabel("Lithium concentration [mol / m^3]")
% grid on;

f3 = figure;
plot(time, j_neg);
hold on;
plot(time, j_pos, 'r');
title("Lithium flux at the electrodes")
xlabel("Time")
ylabel("Lithium flux [mol / m^2 / s]")
grid on;

%f4 = figure;
%plot(time, potse_neg);
%hold on;
%plot(time, potse_pos, 'r');
%title("Potse at the electrodes")
%xlabel("Time")
%ylabel("Potential [V]")
%grid on;

% f5 = figure;
% plot(time, pots_neg);
% % hold on;
% % plot(time, pots_pos, 'r');
% title("Pots at the electrodes")
% xlabel("Time")
% ylabel("Potential [V]")
% grid on;
 
% f6 = figure;
% plot(time, v);
% title("Voltage")
% xlabel("Time")
% ylabel("Potential [V]")
% grid on;
