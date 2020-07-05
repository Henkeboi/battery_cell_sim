function SOC = calculate_SOC(cs, integrator_state, electrode, const)
    F = const.F;
    if electrode == 'neg'
        eps = const.porosity_solid_neg;
        A = const.A_neg;
        L = const.L_neg;
        max_c = const.solid_max_c_neg;
        x100 = const.x100_neg;
        x0 = const.x0_neg;
    elseif electrode == 'pos'
        eps = const.porosity_solid_pos;
        A = const.A_pos;
        L = const.L_pos;
        max_c = const.solid_max_c_pos;
        x100 = const.x100_pos;
        x0 = const.x0_pos;
    else
        error("Bad electrode selection")
    end
    cs_avg = cs - integrator_state / (eps * A * F * L); % Calculate average electrode concentration.
    SOC = (cs_avg / max_c - x0) / (x100 - x0); % Calculate SOC from average concentration.
end
