function [] = tf_ce(cse, x_coordinates, T_len, sampling_f, const)
    [Rct_n, Rfilm_n, Rse_n, Rs_n, as_n, L_n, F_n, Ds_n, A_n, alpha_n, sigma_n, kappa_n, beta_n, eps_s_n, cs0_n, Uocv_d_n] = get_electrode_constants(0, 'neg', const);
    [Rct_p, Rfilm_p, Rse_p, Rs_p, as_p, L_p, F_p, Ds_p, A_p, alpha_p, sigma_p, kappa_p, beta_p, eps_s_p, cs0_p, Uocv_d_p] = get_electrode_constants(0, 'pos', const);
    L_s = const.L_sep;
    L_tot = L_n + L_s + L_p;
    eps_e_n = const.eps_e_neg;
    eps_e_s = const.eps_e_sep;
    eps_e_p = const.eps_e_pos;
    De_n = const.Deeff_neg;
    De_s = const.Deeff_sep;
    De_p = const.Deeff_pos;
    
    syms x lambda w1 w2 k1 k3 k4 k5 k6;
        
    w1 = L_n * sqrt(lambda * eps_e_n / De_n);
    w2 = L_n * sqrt(lambda * eps_e_s / De_s);
    w3 = (L_n + L_s) * sqrt(lambda * eps_e_s / De_s);
    w4 = (L_n + L_s) * sqrt(lambda * eps_e_p / De_p);

    % Solve in terms of k1
    eq1 = k1 * cos(w1) == k3 * cos(w2) + k4 * sin(w2);
    eq2 = -De_n * k1 * w1 * sin(w1) == De_s * (-k3 * w2 * sin(w2) + k4 * w2 * cos(w2));
 
    % Solve in terms of k3 and k4
    eq3 = k3 * cos(w3) + k4 * sin(w3) == k5 * cos(w4) + k6 * sin(w4);
    eq4 = De_s * (-k3 * w3 * sin(w3) / (L_n + L_s) + k4 * w3 * cos(w3) / (L_n + L_s)) == De_p * (-k5 * w4 * sin(w4) / (L_n + L_s) + k6 * w4 * cos(w4) / (L_n + L_s));
    S_34 = solve([eq3, eq4], [k3, k4]);
    
    phi_n = @(x, lambda, k1) k1 * cos(sqrt(lambda * eps_e_n / De_n) * x) * ((0 <= x) & (x < L_n));
    phi_s = @(x, lambda, k3, k4) (k3 * cos(sqrt(lambda * eps_e_s / De_s) * x) + k4 * sin(sqrt(lambda * eps_e_s / De_s) * x)) * ((L_n <= x) & (x < L_n + L_s));
    phi_p = @(x, lambda, k5, k6) (k5 * cos(sqrt(lambda * eps_e_p / De_p) * x) + k6 * sin(sqrt(lambda * eps_e_p / De_p) * x)) * ((L_n + L_s <= x) & (x <= L_tot));
    phi = @(x, lambda, k1, k3, k4, k5, k6) phi_n(x, lambda, k1) + phi_s(x, lambda, k3, k4) + phi_p(x, lambda, k5, k6);

    fun_eps_e_n = @(x) eps_e_n * ((0 <= x) & (x < L_n));
    fun_eps_e_s = @(x) eps_e_s * ((L_n <= x) & (x < L_n + L_s));
    fun_eps_e_p = @(x) eps_e_p * ((L_n + L_s <= x) & (x <= L_tot));
    fun_eps_e = @(x) fun_eps_e_n(x) + fun_eps_e_s(x) + fun_eps_e_p(x)

end
