sampling_freq = 200;
T = 1 / sampling_freq;
min_T_len = 50;
num_samples = 2 ^ (ceil(log2(min_T_len * sampling_freq)));
freq_vector = 1 : num_samples - 1;
s = (2j / T) * tan(pi * freq_vector / num_samples);

num = s / diffusivity_pos * radius_pos ^ 2.0 + 3.0 - 3.0 * radius_pos * sqrt(s / diffusivity_pos) .* coth(radius_pos * sqrt(s / diffusivity_pos));
den = s .* radius_pos .* (1.0 - radius_pos .* sqrt(s / diffusivity_pos) .* coth(radius_pos .* sqrt(s / radius_pos)));
H_d = num ./ den; % TODO: Add J(z, s)
h_d = real(ifft(H_d)) * sampling_freq;
hstep = T * cumsum(h_d);

step_size = 0.01;
time_d = 1 : step_size : min_T_len;
td = T * (1 : num_samples - 1);
pulse = [0 diff(interp1(td, hstep, time_d))];

hankel_size = 1;
hankel_matrix = hankel(pulse(:, 2 : end));
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

A_estimate = [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)];
B_estimate = [extended_controllability(:, 1); 0.1];
C_estimate = [extended_observability(1, :), 0];
D_estimate = [pulse(1)];
% stem(time_d, pulse);

debiased_surface_concentration = zeros(hankel_size + 1, 1);
u = 1 * ones(hankel_size + 1, 1);
for i = 0 : 1000
    next_debiased_surface_concentration = A_estimate .* debiased_surface_concentration + B_estimate .* u; 
    y = C_estimate * debiased_surface_concentration + D_estimate * u;
    debiased_surface_concentration = next_debiased_surface_concentration;
end



