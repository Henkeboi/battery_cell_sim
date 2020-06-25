function [] = dra(transfer_function, res0, sampling_freq, T_len, const)
    % sampling_freq and T_len needs to be the same as in the transfer function.
    my_TF = transfer_function(:, 1); % Temp Simpify
    
    T = 1 / sampling_freq;
    num_samples = 2 ^ (ceil(log2(T_len * sampling_freq)));

    my_tf  = real(ifft(my_TF)) * sampling_freq; 
    h_step = T * cumsum(my_tf); 
    td = T * (0 : num_samples - 1);

    T_shifted = 0.01;
    time_vector = 0 : T_shifted : T_len;

    h_pulse = [0 diff(interp1(td, h_step, time_vector))];

    hankel_size = 3;
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

    A_est = [pinv(extended_observability) * H_shifted * pinv(extended_controllability) zeros(system_order, 1); 1 zeros(1, system_order)];
    B_est = [extended_controllability(:, 1); T_shifted]; 
    C_est = [extended_observability(1, :), res0(1)]; % z(1) choosen
    D_est = [0];

    sys = ss(A_est, B_est, C_est, D_est, T_shifted);
    gain = dcgain(sys);

    sys_scaled = ss(A_est, B_est / gain, C_est, D_est, T_shifted);

    steps = 1000;
    li_flux_neg_vector = zeros(steps, 1);
    X_li_flux_neg = zeros(size(A_est, 2), 1);
    for current_step = 1 : steps
        X_li_flux_neg = A_est * X_li_flux_neg + B_est ;
        Y_li_flux_neg = C_est * X_li_flux_neg + D_est ;
        if current_step == 1
            flux_neg_vector(current_step, 1) = 0;
        else
            li_flux_neg_vector(current_step, 1) = li_flux_neg_vector(current_step - 1, 1) + Y_li_flux_neg;
        end
    end
    t_vector = 1 : size(li_flux_neg_vector, 1);
    plot(t_vector, li_flux_neg_vector)

    fprintf("\n\n");


    % tf_pulse = real(ifft(my_tf)) * sampling_freq;
    % T = 0.5;
    % N = 2 ^ (ceil(log2(sampling_freq * T_len)));
    % td = T * (0 : N - 1); 
    % T_sample = 1;
    % t_d = 0 : T_sample : T_len;
    % tf_step = T * cumsum(tf_pulse);
    % tf_d = [0 diff(interp1(td, tf_step, t_d))];
    % hankel_size = 5;
    % big_hankel = hankel(tf_d(1, 2 : end));
    % H = big_hankel(1 : hankel_size, 1 : hankel_size);
    % H_shifted =  big_hankel(2 : hankel_size + 1, 1 : hankel_size);
    % [U S V] = svds(H);
    % sys_order = rank(H);
    % S = S(1 : sys_order, 1 : sys_order);
    % U = U(:, 1 : sys_order);
    % V = V(:, 1 : sys_order);
    % sigma = S .^ 0.5;
    % extended_obs = U * sigma;
    % extended_con = sigma * V';
    % A_estimate = [pinv(extended_obs) * H_shifted * pinv(extended_con) zeros(sys_order, 1); 1 zeros(1, sys_order)];
    % B_estimate = [extended_con(:, 1); T_sample];
    % C_estimate = [extended_obs(1, :), res0];
    % D_estimate = [0];
    % 
    % sys = ss(A_estimate, B_estimate, C_estimate, D_estimate, T_sample);
    % gain = dcgain(sys);
    % sys_scaled = ss(A_estimate, B_estimate / gain, C_estimate, D_estimate, T_sample);


    % sampling_f = 200;
    % T = 1 / sampling_f;
    % min_T_len = 1;
    % num_samples = 2 ^ (ceil(log2(sampling_f * min_T_len)));
    % f_vector = 0 : num_samples - 1;
    % s = zeros(1, size(f_vector, 2));
    % for i = 1 : size(f_vector, 2)
    %     s(1, i) = 2j * sampling_f * tan(pi * f_vector(i) / num_samples);
    % end
    %
    % cse_neg = 2000;
    % Rct = const.R_ct_neg;
    % Rfilm = const.R_film_neg;
    % Rtot = Rct + Rfilm;
    % Rs = const.radius_neg;
    % L = const.L_neg;
    % F = const.F;
    % Ds = const.diffusivity_neg;
    % Acell = const.A_neg;
    % as = const.alpha_neg;
    % sigma_eff = const.sigma_neg;
    % kappa_eff = const.kappa_neg;
    % beta = Rs * sqrt(s / Ds);
    % eps = const.porosity_solid_neg;
    % locs = [0.1 0.3 0.6 1];

    % dudc = calculate_ocv_derivative_neg(cse_neg, const); 
    % beta = Rs * sqrt(s./Ds);          % Jacobsen-West beta
    % % Compute (scaled) transfer function of C_{s,e}(s)/J(s) [Jacobsen-West]  
    % cse_j = Rs/Ds./(1-beta.*coth(beta)); 
    % % Compute unitless impedance ratio term (nu_inf as s->infinity)
    % nu = L * sqrt((as/sigma_eff + as/kappa_eff)./(Rtot + dudc.*cse_j/F));
    % 
    % % Now, we are ready to actually calculate the transfer function
    % len = length(locs); locs = locs(:);
    % j_tf = zeros(len,length(s)); % Initialize output to zero
    % % Initialize "D" variables and transfer-function names variable
    % Dterm = zeros(len,1); Dstr = cell(len,1); names = cell(len,1);
    % for n=1:len,
    %   z = locs(n);
    %   j_tf(n,:) = nu./(as*F*L*Acell*(kappa_eff+sigma_eff).*sinh(nu)).*...
    %               (sigma_eff*cosh(nu.*z) + kappa_eff*cosh(nu.*(z-1)));
    %   if any(isnan(j_tf)),
    %     error('ERROR (tf_j): At least one value computes to NaN');
    %   end
    % end
  
    % % value of tf at s->0. Solved with Mathematica.
    % tf0 = 1/(Acell*as*F*L);
    % j_tf(:, s==0) = tf0(:, ones(size(find(s==0))));


    % tf_pulse = real(ifft(j_tf(1, :))) * sampling_freq;
end
