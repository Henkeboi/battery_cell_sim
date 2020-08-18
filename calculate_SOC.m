function SOC = calculate_SOC(cs, integrator_state, electrode, const)
    if isnan(integrator_state)
        error("Input is NAN")
    end

    % Calculates SOC by subtracting lithium flux out of the particles from initial concentration. 
    F = const.F;
    if electrode == 'neg'
        eps = const.eps_s_neg;
        A = const.A_neg;
        L = const.L_neg;
        max_c = const.solid_max_c_neg;
        x100 = const.x100_neg;
        x0 = const.x0_neg;
    elseif electrode == 'pos'
        eps = const.eps_s_pos;
        A = const.A_pos;
        L = const.L_pos;
        max_c = const.solid_max_c_pos;
        x0 = const.x0_pos;
        x100 = const.x100_pos;
    else
        error("Bad electrode selection")
    end
    cs_avg = cs - integrator_state / (eps * A * F * L); % Calculate average electrode concentration.
    if electrode == 'neg'
        SOC = (cs_avg / max_c - x0) / (x100 - x0); % Calculate SOC from average concentration.
    else
        SOC = 1 + (cs_avg / max_c - x0) / (x100 - x0); % Calculate SOC from average concentration.
    end
end
