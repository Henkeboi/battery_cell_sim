clearvars;
disp("Loading data.")
run('parameters.m');
run('read_load_cycle.m')

%% Calculations for the negative eletrode at z = 0.
cs0_neg = const.solid_max_c_neg * const.x100_neg;
cs0_pos = const.solid_max_c_pos * const.x0_pos;
cse0_neg = cs0_neg / const.radius_neg;
cse0_pos = cs0_neg / const.radius_pos;

disp("Running blended Lithium surface concentration dra.")
cse_neg_sampling_f = 400;
cse_neg_T_len = 5;
cse_neg_blender = Blender(0.01, @tf_cse, [0], 'neg', const);
cse_neg_blender.create_models(cse_neg_T_len, cse_neg_sampling_f);
cse_neg_blender.sort();

cse_pos_sampling_f = 400;
cse_pos_T_len = 5;
cse_pos_blender = Blender(0.01, @tf_cse, [0], 'pos', const);
cse_pos_blender.create_models(cse_pos_T_len, cse_pos_sampling_f);
cse_pos_blender.sort();

disp("Running blended lithium flux dra.")
j_neg_sampling = 400;
j_neg_T_len = 5;
j_neg_blender = Blender(0.01, @tf_j, [0], 'neg', const);
j_neg_blender.create_models(j_neg_T_len, j_neg_sampling);
j_neg_blender.sort();

j_pos_sampling = 400;
j_pos_T_len = 5;
j_pos_blender = Blender(0.01, @tf_j, [0], 'pos', const);
j_pos_blender.create_models(j_pos_T_len, j_pos_sampling);
j_pos_blender.sort();

disp("Running blended potse dra.")
potse_neg_sampling = 400;
potse_neg_T_len = 5;
potse_neg_blender = Blender(0.01, @tf_potse, [0], 'neg', const);
potse_neg_blender.create_models(potse_neg_T_len, potse_neg_sampling);
potse_neg_blender.sort();

potse_pos_sampling = 400;
potse_pos_T_len = 5;
potse_pos_blender = Blender(0.01, @tf_potse, [0], 'pos', const);
potse_pos_blender.create_models(potse_pos_T_len, potse_pos_sampling);
potse_pos_blender.sort();

disp("Simulating.")
z_neg = zeros(size(load_cycle, 1), 1);
z_pos = zeros(size(load_cycle, 1), 1);
cse_neg = zeros(size(load_cycle, 1), 1);
cse_pos = zeros(size(load_cycle, 1), 1);
j_neg = zeros(size(load_cycle, 1), 1);
j_pos = zeros(size(load_cycle, 1), 1);
potse_neg = zeros(size(load_cycle, 1), 1);
potse_pos = zeros(size(load_cycle, 1), 1);
time = zeros(size(load_cycle, 1), 1);
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
    cse_neg(i) = cse0_neg + cse_neg_Y;
    z_neg(i) = calculate_SOC(cs0_neg, cse_neg_X(cse_neg_integrator_index), 'neg', const); 
    SOC_neg = z_neg(i);

    [cse_pos_X, cse_pos_Y, cse_pos_integrator_index] = cse_pos_blender.step(U, SOC_pos);
    cse_pos(i) = cse0_pos + cse_pos_Y;
    z_pos(i) = calculate_SOC(cs0_pos, cse_pos_X(cse_pos_integrator_index), 'pos', const); 
    SOC_pos = z_pos(i);
    
    [j_neg_X, j_neg_Y, j_neg_integrator_index] = j_neg_blender.step(U, SOC_neg);
    j_neg(i) = j_neg_Y;

    [j_pos_X, j_pos_Y, j_pos_integrator_index] = j_pos_blender.step(U, SOC_pos);
    j_pos(i) = j_pos_Y;
 
    [potse_neg_X, potse_neg_Y, potse_neg_integrator_index] = potse_neg_blender.step(U, SOC_neg);
    potse_neg(i) = potse_neg_Y; % TODO: Debieas and debug.

    [potse_pos_X, potse_pos_Y, potse_pos_integrator_index] = potse_pos_blender.step(U, SOC_pos);
    potse_pos(i) = potse_pos_Y; % TODO: Debieas and debug.
    

    time(i) = time_acc;
    time_acc = time(i) + delta_time;

    if isnan(z_neg(i))
        disp("Nan z_neg")
        break;
    end
    if isnan(z_pos(i))
        disp("Nan z_pos")
        break;
    end
    if isnan(cse_neg(i))
        disp("Nan cse_neg")
        break;
    end
    if isnan(cse_pos(i))
        disp(i);
        disp("nan cse_pos")
        break;
    end
    if isnan(j_neg(i))
        disp("nan j_neg")
        break;
    end
    if isnan(j_pos(i))
        disp("nan j_neg")
        break;
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

f2 = figure;
plot(time, cse_neg);
hold on;
plot(time, cse_pos, 'r');
title("Average surface concentration")
xlabel("Time")
ylabel("Lithium concentration [mol / m^3]")
grid on;

f3 = figure;
plot(time, j_neg);
hold on;
plot(time, j_pos, 'r');
title("Lithium flux at the electrodes")
xlabel("Time")
ylabel("Lithium flux [mol / m^2 / s]")
grid on;

f4 = figure;
plot(time, potse_neg);
hold on;
plot(time, potse_pos, 'r');
title("Potse at the electrodes")
xlabel("Time")
ylabel("Potential [V]")
grid on;








