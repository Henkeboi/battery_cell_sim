% Cell parameters
const.R = 8.31446261815324;
const.F = 96485.3329;
const.L_neg = 128e-6;
const.L_pos = 190e-6;
const.L_sep = 76e-6;
const.transfer_number_neg = 0.363;
const.transfer_number_sep = 0.363;
const.transfer_number_pos = 0.363;
const.A_neg = 1;
const.A_pos = 1;
const.A = 1;
const.radius_neg = 12.5e-6;
const.radius_pos = 8.5e-6;
const.eps_s_neg = 0.471;
const.eps_s_pos = 0.297;
const.eps_e_neg = 0.357;
const.eps_e_sep = 0.724;
const.eps_e_pos = 0.444;
const.brug_neg = 1.5;
const.brug_sep = 1.5;
const.brug_pos = 1.5;
const.reaction_rate_neg = 1.94e-11;
const.reaction_rate_pos = 2.16e-11;

const.kappa_norm_neg = 2.2842e-5;
const.kappa_norm_pos = 2.2073e-5;
const.sigma_norm_neg = 100;
const.sigma_norm_pos = 3.8;

const.kappa_eff_neg = const.kappa_norm_neg * (const.eps_s_neg ^ const.brug_neg);
const.kappa_eff_pos = const.kappa_norm_pos * (const.eps_s_pos ^ const.brug_pos);
%const.kappa_eff_neg = const.reaction_rate_neg * (const.eps_s_neg ^ const.brug_neg);
%const.kappa_eff_pos = const.reaction_rate_pos * (const.eps_s_pos ^ const.brug_pos);

const.sigma_eff_neg = const.sigma_norm_neg * (const.eps_s_neg ^ const.brug_neg);
const.sigma_eff_pos = const.sigma_norm_pos * (const.eps_s_pos ^ const.brug_pos);

const.as_neg = 3 * const.eps_s_neg / const.radius_neg;
const.as_pos = 3 * const.eps_s_pos / const.radius_pos;
const.x100_neg = 0.4;
const.x100_pos = 0.9;
const.x0_neg = 0.1;
const.x0_pos = 0.5;
const.solid_max_c_neg = 10000;
const.solid_max_c_pos = 10000;
const.cs0_neg = const.solid_max_c_neg * const.x100_neg;
const.cs0_pos = const.solid_max_c_pos * const.x0_pos;
const.ce0_neg = 10000;
const.ce0_sep = 10000;
const.ce0_pos = 10000;
const.diffusivity_neg = 3.9e-14;
const.diffusivity_pos = 1e-13;
const.De_neg = 7.5e-11;
const.De_sep = 7.5e-11;
const.De_pos = 7.5e-11;
const.Deeff_neg = const.De_neg * const.eps_e_neg ^ const.brug_neg;
const.Deeff_sep = const.De_sep * const.eps_e_sep ^ const.brug_sep;
const.Deeff_pos = const.De_pos * const.eps_e_pos ^ const.brug_pos;
const.alpha_neg = 0.5;
const.alpha_pos = 0.5;
const.temp = 25;

const.j0_neg = const.reaction_rate_neg * const.ce0_neg ^ (1 - const.alpha_neg) * (const.solid_max_c_neg - const.cs0_neg) ^ (1 - const.alpha_neg) * const.cs0_neg ^ const.alpha_neg; 
const.j0_pos = const.reaction_rate_pos * const.ce0_pos ^ (1 - const.alpha_pos) * (const.solid_max_c_pos - const.cs0_pos) ^ (1 - const.alpha_pos) * const.cs0_pos ^ const.alpha_pos; 
%const.j0_neg = const.kappa_norm_neg * const.ce0_neg ^ (1 - const.alpha_neg) * (const.solid_max_c_neg - const.cs0_neg) ^ (1 - const.alpha_neg) * const.cs0_neg ^ const.alpha_neg; 
%const.j0_pos = const.kappa_norm_pos * const.ce0_pos ^ (1 - const.alpha_pos) * (const.solid_max_c_pos - const.cs0_pos) ^ (1 - const.alpha_pos) * const.cs0_pos ^ const.alpha_pos; 

const.R_film_neg = 0.0;
const.R_film_pos = 0.0;
const.R_ct_neg = const.R * const.temp / (const.j0_neg * const.F ^ 2);
const.R_ct_pos = const.R * const.temp / (const.j0_pos * const.F ^ 2);
const.R_solid_electrolyte_neg = const.R_film_neg + const.R_ct_neg;
const.R_solid_electrolyte_pos = const.R_film_pos + const.R_ct_pos;
