clear
close all
clc

% Impeller
eta_p = 0.746;
P01 = 65;
P02 = 1602;
pho1 = 71.4;
Q = 1.5898;
m_dot = 1806.205;
Dt2 = 0.4953;
b2 = 0.06858;
beta2 = deg2rad(25);
omega = 575.12;

% Inducer
P1 = 448159;
Psat = 1013.25;
pho = 1141;
rapportoD = 0.81;
theta = deg2rad(10);

[H_m_LOX, delta_h_LOX, vr2_LOX, psi_LOX, P_LOX, ds_LOX, ns_LOX, prod_LOX, NPSP_LOX, NPSH_LOX, sigma_LOX, phi_t_LOX, Cm_LOX, Ss_LOX, tao_LOX, Zt_LOX] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, b2, beta2, omega, P1, Psat, pho, rapportoD, theta);

%%
U2 = omega*(Dt2/2);
Dt1 = rapportoD * Dt2;
R1 = Dt1/2;
U1 = omega*R1;

deltah0 = P_LOX / m_dot;

% ipotesi di vt1 = 0

vt2_controllo = (deltah0)/U2;
w2r = (deltah0/U2.^2 - 1) * (U2 / tan(beta2));
w2t = w2r*tan(beta2);
% w2r = abs((deltah0/U2.^2 - 1) * (U2 / tan(beta2)));
% w2t = abs(w2r*tan(beta2));
v2t = (U2 + w2t);
v2r = w2r;

vt1 = (U2*v2t - deltah0) / U1; % ipotesi verificata
wt1 = U1;

figure
grid on;
grid minor;
hold on;
quiver(0, 0, U1, 0);
quiver(U1, 0, v2t, v2r);
quiver(U1, 0, w2t, w2r);
quiver(U1, v2r, U2, w2t);
legend('U1', 'v2', 'w2', 'U2');

