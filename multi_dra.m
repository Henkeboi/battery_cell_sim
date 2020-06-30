function [A, B, C, D, integrator_index] = multi_dra(A_est, B_est, C_est, D_est, Ts, res0)
    A = [A_est zeros(size(A_est, 1), 1); zeros(1, size(A_est, 1)) 1];
    B = [B_est; Ts];
    C = [C_est res0];
    D = [0];
    integrator_index = size(A, 2);
end
