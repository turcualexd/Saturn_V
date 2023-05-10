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

[H_m_LOX, delta_h_LOX, psi_LOX, P_LOX, ds_LOX, ns_LOX, prod_LOX, NPSP_LOX, NPSH_LOX, sigma_LOX, phi_t_LOX, Cm_LOX, Ss_LOX, tao_LOX, Zt_LOX] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, b2, beta2, omega, P1, Psat, pho, rapportoD, theta);

U2 = omega*(Dt2/2);
Dt1 = rapportoD * Dt2;
R1 = Dt1/2;
U1 = omega*R1;

%%
g = 32.2;
% Q 4070 lmb/sec
Q = 25200;

% H_m_LOX_f
H = (144.*(P02 - P01) ./ pho1);
H_e = 0.30*H;
H_imp = H - H_e;

Dt2_f = Dt2/0.3048;
R2 = Dt2_f/2;
Dt1_f = rapportoD * Dt2_f;
R1 = Dt1_f/2;

A1 = pi*(R1.^2)/2;
A2 = pi*(R2.^2)/2;

cm1 = Q / (448.8*A1);
cm2 = Q / (448.8*A2);

rpm = 5488;

u1 = U1 * 3.281;
u2 = U2 * 3.281;

cu2 = u2 - (cm2 / tan(beta2));
cu1 = cu2/0.75;

c1 = sqrt(cu1.^2 + cm1.^2);
c2 = sqrt(cu2.^2 + cm2.^2);

vu1 = sqrt((u1 - cu1).^2 + cm1.^2);
vm1 = cm1;
v1 = sqrt(vu1.^2 + vm1.^2);

vu2 = sqrt((u2 - cu2).^2 + cm2.^2);
vm2 = cm2;
v2 = sqrt(vu2.^2 + vm2.^2);

u1_m = u1/3.281;
u2_m = u2/3.281;

cm1_m = cm1/3.281;
cm2_m = cm2/3.281;
cu2_m = cu2/3.281;
cu1_m = cu1/3.281;
c1_m = c1/3.281;
c2_m = c2/3.281;

vu1_m = vu1/3.281;
vm1_m = vm1/3.281;
v1_m = v1/3.281;
vu2_m = vu2/3.281;
vm2_m = vm2/3.281;
v2_m = v2/3.281;

figure
grid on;
grid minor;
hold on;
q1 = quiver(0, 0, u1_m, 0, 1);
q1.MaxHeadSize = 0.08;
q2 = quiver(u1_m, 0, -vu1_m, vm1_m, 1);
q2.MaxHeadSize = 0.1;
q3 = quiver(0, 0, u1_m-vu1_m, cm1_m, 1);
q3.MaxHeadSize = 0.1;

q4 = quiver(u1_m, 0, u2_m, 0, 1);
q4.MaxHeadSize = 0.08;
q5 = quiver(u1_m+u2_m, 0, -vu2_m, vm2_m, 1);
q5.MaxHeadSize = 0.1;
q6 = quiver(u1_m, 0, u2_m-vu2_m, cm2_m, 1);
q6.MaxHeadSize = 0.1;

legend('u1', 'v1', 'c1', 'u2', 'v2', 'c2');
title('Triangolo di velocità POMPA LOX');
