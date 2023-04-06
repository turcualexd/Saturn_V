function [rp1_tank, lox_tank] = saturn_v_tank_sizes(rp1_mass, lox_mass, rp1_density, lox_density, stage_diameter, rp1_pressure, lox_pressure)
% il raggio, la lunghezza e lo spessore delle pareti.
% calcola i volumi dei serbatoi
rp1_volume = rp1_mass / rp1_density;
lox_volume = lox_mass / lox_density;

% calcola le dimensioni dei serbatoi
rp1_radius = stage_diameter / 2;
rp1_length = rp1_volume / (pi * rp1_radius^2);
rp1_thickness = rp1_pressure * rp1_radius / (2 * 65000);
rp1_tank = [rp1_radius, rp1_length, rp1_thickness];

lox_radius = stage_diameter / 2;
lox_length = lox_volume / (pi * lox_radius^2);
lox_thickness = lox_pressure * lox_radius / (2 * 65000);
lox_tank = [lox_radius, lox_length, lox_thickness];
end
