clear
close all
clc

%% LOX

% Impeller
eta_p = 0.746;
P01 = 65;
P02 = 1602;
pho1 = 71.4;
Q = 1.5898;
m_dot = 1806.205;
Dt2 = 0.4953;
b2 = 0.06858;
beta2 = deg2rad(35);
omega = 575.12;

% Inducer
P1 = 448159;
Psat = 1013.25;
pho = 1141;
rapportoD = 0.81;
theta = deg2rad(10);

[H_m_LOX, delta_h_LOX, psi_LOX, P_LOX, ds_LOX, ns_LOX, prod_LOX, NPSP_LOX, NPSH_LOX, sigma_LOX, phi_t_LOX, Cm_LOX, Ss_LOX, tao_LOX, Zt_LOX] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, b2, beta2, omega, P1, Psat, pho, rapportoD, theta);

%% RP-1

% Impeller
eta_p = 0.760;
P01 = 45;
P02 = 1860;
pho1 = 50.5;
Q = 0.9833;
m_dot = 796.5082;
Dt2 = 0.59436;
b2 = 0.04318;
beta2 = deg2rad(35);
omega = 575.12;

% Inducer
P1 = 310264;
Psat = 137895;
pho = 810;
rapportoD = 0.81;
theta = deg2rad(10);

[H_m_RP1, delta_h_RP1, psi_RP1, P_RP1, ds_RP1, ns_RP1, prod_RP1, ...
    NPSP_RP1, NPSH_RP1, sigma_RP1, phi_t_RP1, Cm_RP1, Ss_RP1, tao_RP1, Zt_RP1] = ...
    turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, b2, beta2, omega, P1, Psat, pho, rapportoD, theta);