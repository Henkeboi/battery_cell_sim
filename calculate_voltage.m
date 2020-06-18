% x-coordiante = 0 at neg electrode current collector and L at pos electrode current collector.
% z-coordinate = 0 at the current collector and z = 1 at the electrolyte interface.
% U in z coordinates eg z = x / L.
% voltage(t) = Ø_pos(0, t) - Ø_neg(0, t) in z coordinates.
% We need to know V to debias debiased_Ø_neg.
% Linearization of the Butler-Volmer equiation: 

% TODO: Currently using the initial values for concentration and j.
n_pos = 2 * R * T / F * asinh(j0_pos / (2 * reaction_rate_constant_pos * sqrt(initial_c_pos_electrolyte * (max_c_pos - initial_c_pos_electrolyte) * initial_c_pos_electrolyte)));
n_neg = 2 * R * T / F * asinh(j0_neg / (2 * reaction_rate_constant_neg * sqrt(initial_c_neg_electrolyte * (max_c_neg - initial_c_neg_electrolyte) * initial_c_neg_electrolyte)));
voltage = F * (R_film_pos * j0_pos - R_film_neg * j0_neg); % TODO finish expression


