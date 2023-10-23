clear
close all
clc

% Impeller
eta_p = 0.760;
P01 = 45;
P02 = 1860;
pho1 = 50.5;
Q = 0.9617;
m_dot = 796.5082;
Dt2 = 0.59436;
b2 = 0.04318;
beta2 = deg2rad(25);
omega = 574.020;

% Inducer
P1 = 310264;
Psat = 137895;
pho = 808.9324;
rapportoD = 0.81;
theta = deg2rad(10);

[H_m_RP1, delta_h_RP1, psi_RP1, P_RP1, ds_RP1, ns_RP1, prod_RP1, NPSP_RP1, NPSH_RP1, sigma_RP1, phi_t_RP1, Cm_RP1, Ss_RP1, tao_RP1, Zt_RP1] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, b2, beta2, omega, P1, Psat, pho, rapportoD, theta);

%%
U2 = omega*(Dt2/2);
Dt1 = rapportoD * Dt2;
R1 = Dt1/2;
U1 = omega*R1;

deltah0 = P_RP1 / m_dot;

% ipotesi di vt1 = 0

vt2_controllo = (deltah0)/U2;
w2r = (deltah0/U2.^2 - 1) * (U2 / tan(beta2));
w2t = w2r*tan(beta2);
% w2r = abs((deltah0/U2.^2 - 1) * (U2 / tan(beta2)));
% w2t = abs(w2r*tan(beta2));
v2t = (U2 + w2t);
v2r = w2r;

vt1 = (U2*v2t - deltah0) / U1;
wt1 = U1;

figure
grid on;
grid minor;
hold on;
q1 = quiver(0, 0, U1, 0, 1);
q2 = quiver(U1, 0, v2t, v2r, 1);
q3 = quiver(U1, 0, w2t, w2r, 1);
q4 = quiver(w2t+U1, w2r, U2, 0, 1);
legend('U1', 'v2', 'w2', 'U2');
title('Triangolo di velocità POMPA RP-1');
