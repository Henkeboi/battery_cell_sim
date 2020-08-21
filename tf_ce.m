function [] = tf_ce(cse, x_coordinates, T_len, sampling_f, const)
    %eq6 = int(fun_phi_s(x, lambda, k3, k4) .^ 2 * eps_e_s, x, L_n, L_n + L_s) == 1;
    %solve(eq6, k1)
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
    [k3, k4, k5, k6] = get_constants(k1, const);
    digits(30)
    k1 = vpa(k1);
    k3 = vpa(k3);
    k4 = vpa(k4);
    k5 = vpa(k5);
    k6 = vpa(k6);
    lambda = vpa(lambda);

    
    fun_Deeef = @(x) De_n .* ((0 <= x) & (x < L_n)) + De_s.* ((L_n <= x) & (x < L_n + L_s)) + De_p .* ((L_n + L_s <= x) & (x <= L_tot));
    fun_eps = @(x) eps_e_n .* ((0 <= x) & (x < L_n)) + eps_e_s .* ((L_n <= x) & (x < L_n + L_s)) + eps_e_p .* ((L_n + L_s <= x) & (x <= L_tot));
    fun_phi_n = @(x, lambda, k1) k1 .* cos(sqrt(lambda * fun_eps(x) / De_n) .* x);
    fun_phi_s = @(x, lambda, k3, k4) (k3 .* cos(sqrt(lambda * fun_eps(x) / De_s) .* x) + k4 .* sin(sqrt(lambda * fun_eps(x) / De_s) .* x));
    fun_phi_p = @(x, lambda, k5, k6) (k5 .* cos(sqrt(lambda * fun_eps(x) / De_p) .* x) + k6 .* sin(sqrt(lambda * fun_eps(x) / De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot));

    d_phi = diff(fun_phi_n(x, lambda, k1) .* ((0 <= x) & (x < L_n)) + fun_phi_s(x, lambda, k3, k4) .* ((L_n <= x) & (x < L_n + L_s)) + fun_phi_p(x, lambda, k5, k6) .* ((L_n + L_s <= x) & (x <= L_tot)), x);
    
    %A = [fun_Deeef(0) fun_Deeef(L_n) fun_Deeef(L_n + L_s);
    %    fun_Deeef(0) fun_Deeef(L_n) fun_Deeef(L_n + L_s);
    %    fun_Deeef(0) fun_Deeef(L_n) fun_Deeef(L_n + L_s);];
    %V = orth(A);
    %R = V' * A * V;
    %eig(R)

    %return

    tic
    disp("Calculating PHI neg.")
    PHI_n = int(fun_phi_n(x, lambda, k1) .^ 2 * fun_eps, x, 0, L_n);
    disp("Calculating PHI sep.")
    PHI_s = int(fun_phi_s(x, lambda, k3, k4) .^ 2 * fun_eps, x, L_n, L_n + L_s);
    disp("Calculate PHI pos.")
    PHI_p = int(fun_phi_p(x, lambda, k5, k6) .^ 2 * fun_eps, x, L_n + L_s, L_tot);
    eq5 = PHI_n + PHI_s + PHI_p == 1;
    %disp("Solving SL ")
    %solve([eq5], [k1])
    disp("test")
    eq6 = diff(fun_phi_n(x, lambda, k1) .* ((0 <= x) & (x < L_n)) + fun_phi_s(x, lambda, k3, k4) .* ((L_n <= x) & (x < L_n + L_s)) + fun_phi_p(x, lambda, k5, k6) .* ((L_n + L_s <= x) & (x <= L_tot)), x) == 0;
    assume(0 <= lambda & lambda <= 1 & 'real');
    assume(0 <= x & x <= L_tot & 'real')
    assume(0 <= k1 & k1 <= 1 & 'real')
    solve([eq5, eq6], [k1, lambda]) 
    toc

end

function [k3, k4, k5, k6] = get_constants(k1, const)
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
    
    syms x lambda w1 w2 k3 k4 k5 k6;
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
end
