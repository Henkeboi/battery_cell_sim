% Cell parameters
const.R = 8.31446261815324;
const.F = 96485.3329;
const.L_neg = 128e-6;
const.L_pos = 190e-6;
const.L_sep = 76e-6;
const.transfer_neg = 0.363;
const.transfer_sep = 0.363;
const.transfer_pos = 0.363;
const.A_neg = 1;
const.A_pos = 1;
const.A = 1;
const.radius_neg = 12.5e-6;
const.radius_pos = 8.5e-6;
const.resistivity_neg = 100;
const.resistivity_pos = 3.8;
const.porosity_solid_neg = 0.471;
const.porosity_solid_pos = 0.297;
const.porosity_electrolyte_neg = 0.357;
const.porosity_electrolyte_sep = 0.724;
const.porosity_electrolyte_pos = 0.444;
const.brug_neg = 1.5;
const.brug_sep = 1.5;
const.brug_pos = 1.5;
const.reaction_rate_neg = 1.94e-11;
const.reaction_rate_pos = 2.16e-11;
const.kappa_norm_neg = 2.2842e-5; % Normalized reactoin rate.
const.kappa_norm_pos = 2.2073e-5; % Normalized reactoin rate.
const.k_neg =  1.94 * 10e-11; % Reaction_rate_constant.
const.k_pos =  2.16 * 10e-11; % Reaction_rate_constant.
const.kappa_neg = const.reaction_rate_neg * (const.porosity_electrolyte_neg ^ const.brug_neg);
const.kappa_pos = const.reaction_rate_pos * (const.porosity_electrolyte_pos ^ const.brug_pos);
const.sigma_neg = const.resistivity_neg * (const.porosity_solid_neg ^ const.brug_neg);
const.sigma_pos = const.resistivity_pos * (const.porosity_solid_pos ^ const.brug_pos);
const.x100_neg = 0.8;
const.x100_pos = 0.8;
const.x0_neg = 0.2;
const.x0_pos = 0.2;
const.solid_max_c_neg = 26390;
const.solid_max_c_pos = 22860;
const.initial_electrolyte_concentration_neg = 4000;
const.initial_electrolyte_concentration_sep = 2000;
const.initial_electrolyte_concentration_pos = 2000;
const.ce0 = 2000;
const.cs0_neg = 12000;
const.cs0_pos = 12000;
const.diffusivity_neg = 3.9e-14;
const.diffusivity_pos = 1e-13;
const.diffusivity_neg_electorlyte = 7.5e-11;
const.diffusivity_sep = 7.5e-11;
const.diffusivity_pos_electrolyte = 7.5e-11;
const.alpha_neg = 0.5;
const.alpha_pos = 0.5;
const.ionic_conductivity_neg = 1.94e-11;
const.ionic_conductivity_pos = 2.16e-11;
const.temp = 25;
const.j0_neg = const.reaction_rate_neg * const.initial_electrolyte_concentration_neg ^ (1 - const.alpha_neg) * (const.solid_max_c_neg - const.initial_electrolyte_concentration_neg) ^ (1 - const.alpha_neg) * const.initial_electrolyte_concentration_neg ^ const.alpha_neg; 
const.j0_pos = const.reaction_rate_pos * const.initial_electrolyte_concentration_pos ^ (1 - const.alpha_pos) * (const.solid_max_c_pos - const.initial_electrolyte_concentration_pos) ^ (1 - const.alpha_pos) * const.initial_electrolyte_concentration_pos ^ const.alpha_pos; 
const.R_film_neg = 0.0;
const.R_film_pos = 0.0;
const.R_ct_neg = const.R * const.temp / (const.j0_neg * const.F ^ 2);
const.R_ct_pos = const.R * const.temp / (const.j0_pos * const.F ^ 2);
const.R_solid_electrolyte_neg = const.R_film_neg + const.R_ct_neg;
const.R_solid_electrolyte_pos = const.R_film_pos + const.R_ct_pos;
const.R_neg = 12.5e-6;
const.R_pos = 8.5e-6;

% Simulation parameters
const.step_size = 0.1;


% Temporary params
debiased_potential_neg = 0; % Init to 0
initial_potential_pos = 3.8; % TODO: Should be voltage from



