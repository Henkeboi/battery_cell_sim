function [tf_ce, res0, D] = tf_ce(cse, x_location, T_len, freq, M, const)
    L_n = const.L_neg;
    L_p = const.L_pos;
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

    %[k3, k4, k5, k6, k5_simple, k6_simple, eq5] = calculate_constants(k1, const); 
    load('data/constants.mat');
    %[lambda_list] = calculate_lambdas(k3, k4, k5, k6, k5_simple, k6_simple, M, const);
    load('data/lambda_list.mat');
    %[k1_list, k3_list, k4_list, k5_list, k6_list] = eval_constants(lambda_list, eq5, lambda, k1, k3, k4, k5, k6);
    load('data/constants_eval.mat');
    
    tf_ce = 0;
    fun_eps = @(x) eps_e_n .* ((0 <= x) & (x < L_n)) + eps_e_s .* ((L_n <= x) & (x < L_n + L_s)) + eps_e_p .* ((L_n + L_s <= x) & (x <= L_tot));
    for i = 1 : size(lambda_list, 2)
        lambda = lambda_list(1, i);
        k1 = k1_list(1, i);
        k3 = k3_list(1, i);
        k4 = k4_list(1, i);
        k5 = k5_list(1, i);
        k6 = k6_list(1, i);

        [s, j_n_neg] = calculate_jn(cse, freq, T_len, lambda, k1, k5, k6, 'neg', const);
        [s, j_n_pos] = calculate_jn(cse, freq, T_len, lambda, k1, k5, k6, 'pos', const);
        tf_ce_n = (j_n_neg + j_n_pos) ./ (s + lambda);

        fun_phi_neg = @(x) k1 .* cos(sqrt(lambda .* fun_eps(x) ./ De_n) .* x) .* ((0 <= x) & (x < L_n));
        fun_phi_sep = @(x) ((k3 .* cos(sqrt(lambda .* fun_eps(x) ./ De_s) .* x) + k4 .* sin(sqrt(lambda .* fun_eps(x) ./ De_s) .* x))) .* ((L_n <= x) & (x < L_n + L_s)); 
        fun_phi_pos = @(x) ((k5 .* cos(sqrt(lambda .* fun_eps(x) ./ De_p) .* x) + k6 .* sin(sqrt(lambda .* fun_eps(x) ./ De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot)));
        phi_n = @(x) fun_phi_neg(x) + fun_phi_sep(x) + fun_phi_pos(x);
        tf_ce = tf_ce + tf_ce_n  * phi_n(x_location);
    end
    
    res0 = 0;
    D = NaN;

    if any(isnan(tf_ce))
        error("NAN in tf_ce");
    end
end

function [k3, k4, k5, k6, k5_simple, k6_simple, eq5] = calculate_constants(k1, const)
    L_n = const.L_neg;
    L_p = const.L_pos;
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
    eq3 = k5 * cos(w4) + k6 * sin(w4) == k3 * cos(w3) + k4 * sin(w3);
    eq4 = De_p * (-k5 * w4 * sin(w4) + k6 * w4 * cos(w4)) == De_s * (-k3 * w3 *sin(w3) + k4 * w3 * cos(w3));
    S_56 = solve([eq3, eq4], [k5, k6]);
    k5 = S_56.k5;
    k6 = S_56.k6;
    
    fun_eps = @(x) eps_e_n .* ((0 <= x) & (x < L_n)) + eps_e_s .* ((L_n <= x) & (x < L_n + L_s)) + eps_e_p .* ((L_n + L_s <= x) & (x <= L_tot));
    fun_phi_n = @(x, lambda, k1) k1 .* cos(sqrt(lambda .* fun_eps(x) ./ De_n) .* x) .* ((0 <= x) & (x < L_n));
    fun_phi_s = @(x, lambda, k3, k4) ((k3 .* cos(sqrt(lambda .* fun_eps(x) ./ De_s) .* x) + k4 .* sin(sqrt(lambda .* fun_eps(x) ./ De_s) .* x))) .* ((L_n <= x) & (x < L_n + L_s)); 
    fun_phi_p = @(x, lambda, k5, k6) ((k5 .* cos(sqrt(lambda .* fun_eps(x) ./ De_p) .* x) + k6 .* sin(sqrt(lambda .* fun_eps(x) ./ De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot)));

    disp("Calculating PHI neg.")
    step_size = L_n / 50;
    xn_vector = 0 : step_size : L_n;
    PHI_n = @(x, lambda, k1) trapz(xn_vector, fun_phi_n(xn_vector, lambda, k1) .^ 2 .* fun_eps(xn_vector) .* (0 <= x & x <= L_n));
    disp("Calculating PHI sep.")
    step_size = L_s / 50;
    xs_vector = L_n : step_size : L_n + L_s;
    PHI_s = @(x, lambda, k3, k4) trapz(xs_vector, fun_phi_s(xs_vector, lambda, k3, k4) .^ 2 .* fun_eps(xs_vector) .* (L_n <= x & x <= L_n + L_s));
    disp("Calculate PHI pos.")
    step_size = L_p / 50;
    xp_vector = L_n + L_s : step_size : L_tot;
    PHI_p = @(x, lambda, k5, k6) trapz(xp_vector, fun_phi_p(xp_vector, lambda, k5, k6) .^ 2 .* fun_eps(xp_vector) .* (L_n + L_s <= x <= L_tot)); 

    disp("Isolate k1.")
    eq5 = PHI_n(xn_vector, lambda, k1) + PHI_s(xs_vector, lambda, k3, k4) + PHI_p(xp_vector, lambda, k5, k6) == 1;
    eq5 = isolate(eq5, k1);
    
    disp("Substitute k1.")
    k3 = subs(k3, k1, rhs(eq5));
    k4 = subs(k4, k1, rhs(eq5));
    k5 = subs(k5, k1, rhs(eq5));
    k6 = subs(k6, k1, rhs(eq5));
    k5_simple = simplify(subs(k5, k1, rhs(eq5)));
    k6_simple = simplify(subs(k6, k1, rhs(eq5)));

    disp("Save variables.")
    save('data/constants.mat', 'k3', 'k4', 'k5', 'k6', 'k5_simple', 'k6_simple', 'eq5');
end

function [lambda_list] = calculate_lambdas(k3, k4, k5, k6, k5_simple, k6_simple, M, const)
    disp("Calculate lambdas.");
    L_n = const.L_neg;
    L_p = const.L_pos;
    L_s = const.L_sep;
    L_tot = L_n + L_s + L_p;
    eps_e_n = const.eps_e_neg;
    eps_e_s = const.eps_e_sep;
    eps_e_p = const.eps_e_pos;
    De_n = const.Deeff_neg;
    De_s = const.Deeff_sep;
    De_p = const.Deeff_pos;

    syms x lambda k1

    fun_eps = @(x) eps_e_n .* ((0 <= x) & (x < L_n)) + eps_e_s .* ((L_n <= x) & (x < L_n + L_s)) + eps_e_p .* ((L_n + L_s <= x) & (x <= L_tot));
    phi = @(x, L) eval(subs(rhs(eq5), lambda, L)) .* cos(sqrt(L .* fun_eps(x) ./ De_n) .* x) .* ((0 <= x) & (x < L_n)) ...
        + (eval(subs(k3, lambda, L)) .* cos(sqrt(L .* fun_eps(x) ./ De_s) .* x) ...
        + eval(subs(k4, lambda, L)) .* sin(sqrt(L .* fun_eps(x) ./ De_s) .* x)) .* ((L_n <= x) & (x < L_n + L_s)) ...
        + (eval(subs(k5, lambda, L)) .* cos(sqrt(L .* fun_eps(x) ./ De_p) .* x) ...
        + eval(subs(k6, lambda, L)) .* sin(sqrt(L .* fun_eps(x) ./ De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot));
    phi_p = @(x, L) (eval(subs(k5_simple, lambda, L)) .* cos(sqrt(L .* fun_eps(x) ./ De_p) .* x) ...
        + eval(subs(k6_simple, lambda, L)) .* sin(sqrt(L .* fun_eps(x) ./ De_p) .* x)) .* ((L_n + L_s <= x) & (x <= L_tot));

    lambda_list = [];
    h = L_p / 10e5;
    l_b = 0.01;
    i = 1;
    while i < M
        l_a = l_b;
        step_a = 0.01;
        a = (phi_p(L_tot - h, l_a) - phi_p(L_tot, l_a)) / h;
        while a < 0
            a = (phi_p(L_tot - h, l_a + step_a) - phi_p(L_tot, l_a + step_a)) / h;
            l_a = l_a + step_a;
        end

        step_b = 0.1;
        l_b = l_a + step_b;
        b = (phi_p(L_tot - h, l_b) - phi_p(L_tot, l_b)) / h;
        while b > 0
            b = (phi_p(L_tot - h, l_b +  step_b) - phi_p(L_tot, l_b + step_b)) / h;
            l_b = l_b + step_b;
        end 

        l_c = l_a + (l_b - l_a) / 2;
        c = phi_p(L_tot - h, l_c) / h;
        last_l_c = 0;
        skip = false;
        tol = 10e-4;
        while abs(c) > tol && ~skip
            disp(c)
            if c > 0
                a = c;
                l_a = l_c;
            else
                b = c;
                l_b = l_c;
            end
            last_l_c = l_c;
            l_c = l_a + (l_b - l_a) / 2;
            c = phi_p(L_tot - h, l_c) / h;
            if l_c == last_l_c
                step_a = step_a * 10;
                skip = true;
            end
        end
        if ~skip
            lambda_list = [lambda_list l_c];
            i = i + 1;
        end
        l_b = l_b + step_a;
    end
    save('data/lambda_list.mat', 'lambda_list')
end

function [k1_list, k3_list, k4_list, k5_list, k6_list] = eval_constants(lambda_list, eq5, lambda, k1, k3, k4, k5, k6) 
    k1_list = [];
    k3_list = [];
    k4_list = [];
    k5_list = [];
    k6_list = [];

    for i = 1 : size(lambda_list, 2)
        disp(i)
        k1_list = [k1_list eval(subs(rhs(eq5), lambda, lambda_list(i)))];
        k3_list = [k3_list eval(subs(k3, lambda, lambda_list(i)))];
        k4_list = [k4_list eval(subs(k4, lambda, lambda_list(i)))];
        k5_list = [k5_list eval(subs(k5, lambda, lambda_list(i)))]; 
        k6_list = [k6_list eval(subs(k6, lambda, lambda_list(i)))]; 
        %sl_criteria = PHI_n(xn_vector, lambda_list(i), k1_list(i)) + PHI_s(xs_vector, lambda_list(i), k3_list(i), k4_list(i)) + PHI_p(xp_vector, lambda_list(i), k5_list(i), k6_list(i));
        %if sl_criteria > 1.01 or sl_criteria < 0.99
        %    disp("SL_criteria not satisfied.")
        %end
    end
    save('data/constants_eval.mat', 'k1_list', 'k3_list', 'k4_list', 'k5_list', 'k6_list');
end

function [s, j_n] = calculate_jn(cse, sampling_f, T_len, lambda_n, k1, k5, k6, electrode, const)
    eps_e_n = const.eps_e_neg;
    eps_e_s = const.eps_e_sep;
    eps_e_p = const.eps_e_pos;
    De_n = const.Deeff_neg;
    De_s = const.Deeff_sep;
    De_p = const.Deeff_pos;
    L_n = const.L_neg;
    L_p = const.L_pos;
    L_s = const.L_sep;
    L_tot = L_n + L_s + L_p;

    [s, nu] = calculate_s_nu(cse, T_len, sampling_f, electrode, const);
    [Rct, Rfilm, Rse, Rs, as, L, F, Ds, A, alpha, sigma, kappa, beta, eps, cs0, Uocv_d, transfer_number] = get_electrode_constants(s, electrode, const);
    
    if electrode == 'neg'
        w_n_neg = L * sqrt(lambda_n * eps_e_n / De_n);
        j_n = k1 * (1 - transfer_number) * w_n_neg * sin(w_n_neg) * (kappa + sigma .* cosh(nu)) .* nu ./ (A * F * (kappa + sigma) * (w_n_neg * w_n_neg + nu .* nu) .* sinh(nu)) + k1 * (1 - transfer_number) * (kappa + sigma * cos(w_n_neg)) * nu .* nu ./ (A * F * (kappa + sigma) * (w_n_neg * w_n_neg + nu .* nu));
        j_n(1, 1) = k1 * (1 - transfer_number) / (A * F * w_n_neg);
    elseif electrode == 'pos'
        w_n_pos = L * sqrt(lambda_n * eps_e_p / De_p);
        w_n_tot = L_tot * sqrt(lambda_n * eps_e_p / De_p);
        w_n_sep = w_n_tot - w_n_pos;
        delta = (1 + exp(-2 * nu)) ./ (1 - exp(-2 *nu));

        j_n = ((k5 * (1 - transfer_number) * w_n_pos * sin(w_n_sep) .* ((kappa ./ sinh(nu))+ sigma .* delta) .* nu) ...
            - (k5 * (1 - transfer_number) * w_n_pos * sin(w_n_tot) .* (sigma ./ sinh(nu) + kappa .* delta) .* nu) ...
            - k5 * (1 - transfer_number) * (sigma * cos(w_n_sep) + kappa * cos(w_n_tot)) .* nu .* nu ...
            + k6 * (1 - transfer_number) * w_n_pos * cos(w_n_tot) * (sigma ./ sinh(nu) + kappa .* delta) .* nu ...
            - k6 * (1 - transfer_number) * w_n_pos * cos(w_n_sep) * (kappa ./ sinh(nu) + sigma .* delta) .* nu ...
            - k6 * (1 - transfer_number) * (sigma * sin(w_n_sep) + kappa * sin(w_n_tot)) .* nu .* nu) ...
            ./ (A * F * (kappa + sigma) .* (w_n_pos * w_n_pos + nu .* nu));
    
        nu0 = nu(1, 1);
        j_n(1, 1) = k5 * (1 - transfer_number) * w_n_pos * sin(w_n_sep) .* (kappa + sigma .* cosh(nu0)) .* nu0 / (A * F * (kappa + sigma) .* (w_n_pos * w_n_pos + nu0 .* nu0)) ...
                    - k5 * (1 - transfer_number) * w_n_pos * sin(w_n_tot) * (sigma + kappa) / (A * F * (kappa + sigma) * w_n_pos * w_n_pos) ...
                    - k5 * (1 - transfer_number) * (sigma * cos(w_n_sep) + kappa * cos(w_n_tot)) .* nu0 .* nu0 ./ (A * F * (kappa + sigma) * (w_n_pos * w_n_pos + nu0 * nu0)) ...
                    + k6 * (1 - transfer_number) * w_n_pos * cos(w_n_tot) * (sigma + kappa) / (A * F * (kappa + sigma) * w_n_pos * w_n_pos) ...
                    - k6 * (1 - transfer_number) * w_n_pos * cos(w_n_sep) * (kappa + sigma) / (A * F * (kappa + sigma) * w_n_pos * w_n_pos) ...
                    - k6 * (1 - transfer_number) * (sigma * sin(w_n_sep) + kappa * sin(w_n_tot)) .* nu0 .* nu0;
    else
        error("Bad electrode selection.")
    end

    if any(isnan(j_n))
        error("NaN in j_n");
    end
end
