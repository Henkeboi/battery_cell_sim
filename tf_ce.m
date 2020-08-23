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

    syms x lambda k1
    syms w1 w2 k3 k4 k5 k6;
    assume(lambda, 'real');
    assume(0 <= x & x <= L_tot & 'real');

    w1 = L_n * sqrt(lambda * eps_e_n / De_n);
    w2 = L_n * sqrt(lambda * eps_e_s / De_s);
    w3 = (L_n + L_s) * sqrt(lambda * eps_e_s / De_s);
    w4 = (L_n + L_s) * sqrt(lambda * eps_e_p / De_p);

    % Solve for k3 and k4 in terms of k1
    eq1 = k1 * cos(w1) == k3 * cos(w2) + k4 * sin(w2);
    eq2 = -De_n * k1 * w1 * sin(w1) == De_s * (-k3 * w2 * sin(w2) + k4 * w2 * cos(w2));
    S_34 = solve([eq1, eq2], [k3, k4]);
    k3 = S_34.k3;
    k4 = S_34.k4;

    % Solve for k5 and k6 in terms of k3 and k4
    eq3 = k3 * cos(w3) + k4 * sin(w3) == k5 * cos(w4) + k6 * sin(w4);
    eq4 = De_s * (-k3 * w3 * sin(w3) / (L_n + L_s) + k4 * w3 * cos(w3) / (L_n + L_s)) == De_p * (-k5 * w4 * sin(w4) / (L_n + L_s) + k6 * w4 * cos(w4) / (L_n + L_s));
    S_56 = solve([eq3, eq4], [k5, k6]);
    k5 = S_56.k5;
    k6 = S_56.k6;
    
    fun_eps = @(x) eps_e_n .* ((0 <= x) & (x < L_n)) + eps_e_s .* ((L_n <= x) & (x < L_n + L_s)) + eps_e_p .* ((L_n + L_s <= x) & (x <= L_tot));
    fun_phi_n = @(x, lambda, k1) k1 .* cos(sqrt(lambda * fun_eps(x) / De_n) .* x) .* ((0 <= x) & (x < L_n));
    fun_phi_s = @(x, lambda, k3, k4) ((k3 .* cos(sqrt(lambda * fun_eps(x) / De_s) .* x) + k4 .* sin(sqrt(lambda * fun_eps(x) / De_s) .* x))) .* ((L_n <= x) & (x < L_n + L_s)); 
    fun_phi_p = @(x, lambda, k5, k6) ((k5 .* cos(sqrt(lambda * fun_eps(x) / De_p) .* x) + k6 .* sin(sqrt(lambda * fun_eps(x) / De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot)));

    step_size = 10e-7;
    disp("Calculating PHI neg.")
    xn_vector = 0 : step_size : L_n;
    PHI_n = @(x, lambda, k1) trapz(fun_phi_n(xn_vector, lambda, k1) .^ 2 .* fun_eps(xn_vector) .* (0 <= x & x <= L_n));
    disp("Calculating PHI sep.")
    xs_vector = L_n : step_size : L_n + L_s;
    PHI_s = @(x, lambda, k3, k4) trapz(fun_phi_s(xs_vector, lambda, k3, k4) .^ 2 .* fun_eps(xs_vector) .* (L_n <= x & x <= L_n + L_s));
    disp("Calculate PHI pos.")
    xp_vector = L_n + L_s : step_size : L_tot;
    PHI_p = @(x, lambda, k5, k6) trapz(fun_phi_p(xp_vector, lambda, k5, k6) .^ 2 .* fun_eps(xp_vector) .* (L_n + L_s <= x <= L_tot)); 
    eq5 = PHI_n(xn_vector, lambda, k1) + PHI_s(xs_vector, lambda, k3, k4) + PHI_p(xp_vector, lambda, k5, k6) == 1;
    eq5 = isolate(eq5, k1);


    k3 = subs(k3, k1, rhs(eq5));
    k4 = subs(k4, k1, rhs(eq5));
    k5 = subs(k5, k1, rhs(eq5));
    k6 = subs(k6, k1, rhs(eq5));

    l = 1;
    k3 = eval(subs(k3, lambda, l));
    k4 = eval(subs(k4, lambda, l));
    k5 = eval(subs(k5, lambda, l));
    k6 = eval(subs(k6, lambda, l));
    k1 = eval(subs(rhs(eq5), lambda, l));

    %eval(subs(PHI_n(xn_vector, lambda, k1) + PHI_s(xs_vector, lambda, k3, k4) + PHI_p(xp_vector, lambda, k5, k6), lambda, l))
end
