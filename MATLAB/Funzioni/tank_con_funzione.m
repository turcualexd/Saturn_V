clear; clc; close all;
m_LOX = 635029.318; %kg 
m_RP1 = 1441516.552; %kg

rho_LOX = 1141; % kg/m^3
rho_RP1 = 810; % kg/m^3

diametro = 10.1;

p_rp1 = 120; %psi
p_lox = 180; %psi
%[lunghezza, volume] = dimensioni_serbatoi(m_LOX, rho_LOX, m_RP1, rho_RP1, diametro)

[lunghezza_rp1, volume_rp1, lunghezza_lox, volume_lox] = dimensioni_serbatoi2(m_LOX, rho_LOX, m_RP1, rho_RP1, diametro, p_rp1, p_lox)

[rp1_tank, lox_tank] = saturn_v_tank_sizes(m_LOX, rho_LOX, m_RP1, rho_RP1, diametro, p_rp1, p_lox)
