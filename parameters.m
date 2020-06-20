R = 8.31446261815324;
F = 96485.3329;
L_neg = 128e-6;
L_pos = 190e-6;
A_neg = 1;
A_pos = 1;
radius_neg = 12.5e-6;
radius_pos = 8.5e-6;
resistivity_neg = 100;
resistivity_pos = 3.8;
porosity_solid_neg = 0.471;
porosity_solid_pos = 0.297;
porosity_electrolyte_neg = 0.357;
porosity_electrolyte_sep = 0.724;
porosity_electrolyte_pos = 0.444;
brug_neg = 1.5;
brug_pos = 1.5;
reaction_rate_constant_neg = 1.94e-11;
reaction_rate_constant_pos = 2.16e-11;
k_eff_neg = reaction_rate_constant_neg * (porosity_electrolyte_neg ^ brug_neg);
k_eff_pos = reaction_rate_constant_pos * (porosity_electrolyte_pos ^ brug_pos);
o_eff_neg = resistivity_neg * (porosity_solid_neg ^ brug_neg);
o_eff_pos = resistivity_pos * (porosity_solid_pos ^ brug_pos);
solid_max_concentration_neg = 26390;
solid_max_concentration_pos = 22860;
initial_electrolyte_average_concentration_neg = 2000; 
initial_electrolyte_average_concentration_sep = 2000;
initial_electrolyte_average_concentration_pos = 2000;
cse_neg = 1000; % Inital negative electrode surface concentration
diffusivity_neg = 3.9e-14;
diffusivity_pos = 1e-13;
diffusivity_neg_electorlyte = 7.5e-11;
diffusivity_sep = 7.5e-11;
diffusivity_pos_electrolyte = 7.5e-11;
asymmetric_charge_transfer_neg = 0.5;
asymmetric_charge_transfer_pos = 0.5;
ionic_conductivity_neg = 1.94e-11;
ionic_conductivity_pos = 2.16e-11;
temp = 25;
j0_neg = reaction_rate_constant_neg * initial_electrolyte_average_concentration_neg ^ (1 - asymmetric_charge_transfer_neg) * (solid_max_concentration_neg - initial_electrolyte_average_concentration_neg) ^ (1 - asymmetric_charge_transfer_neg) * initial_electrolyte_average_concentration_neg ^ asymmetric_charge_transfer_neg; 
j0_pos = reaction_rate_constant_pos * initial_electrolyte_average_concentration_pos ^ (1 - asymmetric_charge_transfer_pos) * (solid_max_concentration_pos - initial_electrolyte_average_concentration_pos) ^ (1 - asymmetric_charge_transfer_pos) * initial_electrolyte_average_concentration_pos ^ asymmetric_charge_transfer_pos; 

R_film_neg = 0.0;
R_film_pos = 0.0;
R_charge_transfer_neg = R * temp / (j0_neg * F ^ 2);
R_charge_transfer_pos = R * temp / (j0_pos * F ^ 2);
R_solid_electrolyte_neg = R_film_neg + R_charge_transfer_neg;
R_solid_electrolyte_pos = R_film_pos + R_charge_transfer_pos;
R_neg = 12.5e-6;
R_pos = 8.5e-6;

% Temporary params
debiased_potential_neg = 0; % Init to 0
initial_potential_pos = 3.8; % TODO: Should be voltage from
