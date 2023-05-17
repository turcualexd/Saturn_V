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
%cu1 = cu2/0.75;
cu1 = (u2*cu2 - H*g) / u1;

c1 = sqrt(cu1.^2 + cm1.^2);
c2 = sqrt(cu2.^2 + cm2.^2);

vu1 = sqrt((u1 - cu1).^2 + cm1.^2);
vm1 = cm1;
v1 = sqrt(vu1.^2 + vm1.^2);

vu2 = sqrt((u2 - cu2).^2 + cm2.^2);
vm2 = cm2;
v2 = sqrt(vu2.^2 + vm2.^2);

alpha1 = atan(cm1/cu1);
alpha1_deg = rad2deg(alpha1);
beta1 = atan(vm1/vu1);
beta1_deg = rad2deg(beta1);

alpha2 = atan(cm2/cu2);
alpha2_deg = rad2deg(alpha2);
beta2_controllo = atan(vm2/vu2);
beta2_deg_controllo = rad2deg(beta2_controllo);

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
axis equal;

pos1 = [0, 0, u1_m, 0];
quiver2(pos1, 'red');
pos2 = [u1_m, 0, -vu1_m, vm1_m];
quiver2(pos2, 'blue');
pos3 = [0, 0, u1_m-vu1_m, cm1_m];
quiver2(pos3, 'yellow');

pos4 = [u1_m, 0, u2_m, 0];
quiver2(pos4, 'red');
pos5 = [u1_m+u2_m, 0, -vu2_m, vm2_m];
quiver2(pos5, 'blue');
pos6 = [u1_m, 0, u2_m-vu2_m, cm2_m];
quiver2(pos6, 'yellow');

xline(u1_m, '--');
str1 = {'INLET'};
text(45, -20, str1);
str2 = {'OUTLET'};
text(180, -20, str2);

%%
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
q5.MaxHeadSize = 0.3;
q6 = quiver(u1_m, 0, u2_m-vu2_m, cm2_m, 1);
q6.MaxHeadSize = 0.1;

xline(u1_m, '--');
str1 = {'INLET'};
text(45, -10, str1);
str2 = {'OUTLET'};
text(180, -10, str2);

legend('u1', 'v1', 'c1', 'u2', 'v2', 'c2', '');
title('Triangolo di velocità IMPELLER POMPA LOX');

%% Inducer
N = 5488;
Q = 25200;
NPSH_c_LOX = 60;

NS = (N*sqrt(Q)) ./ (H^0.75);
% parametro di Thoma
thoma = NPSH_c_LOX / H;
NSS = NS / (thoma)^0.75;

NSS_impeller = 11000; % ipotesi fatta da noi
% NPSH impeller
NPSH_imp = (N*sqrt(Q) ./ NSS_impeller).^(4/3);

H_inducer_ft = NPSH_imp - NPSH_c_LOX;
perc_ind = (H_inducer_ft / H) * 100;

psi_ind = 0.06; % ipotesi

u_tip_mean = sqrt((H_inducer_ft*32.2)/psi_ind);

d_t = (720/(pi*N)) .* u_tip_mean; % inducer mean tip diameter

L_ind = d_t * 0.4;
d_0t = d_t + 2*L_ind/2 * tan(deg2rad(7)); % tip diameter at the inducer inlet
d_1t = d_t - 2*L_ind/2 * tan(deg2rad(7)); % tip diameter at the inducer outlet

rd = 0.3; %ipotesi

d_h = d_t * rd; % mean hub diameter
d_0h = d_h - 2*L_ind/2 * tan(deg2rad(14)); % hub diameter at the inducer inlet
d_1h = d_h + 2*L_ind/2 * tan(deg2rad(14)); % hub diameter at the inducer outlet

Q_inducer = Q * (1 + 0.032 + 0.0175);

cm0_ind = Q_inducer/ (3.12 * (pi/4) * (d_0t.^2 - d_0h.^2));
cu0_ind = 0; % assunzione
cm1_ind = Q_inducer/ (3.12 * (pi/4) * (d_1t.^2 - d_1h.^2));

d0 = sqrt((d_0t.^2 + d_0h.^2) / 2); % inducer mean effective diameter inlet

% inducer peripheral velocity at d0
u0_ind = pi*(7000/720)*d0;

d1 = sqrt((d_1t.^2 + d_1h.^2) / 2); % inducer mean effective diameter outlet

% inducer peripheral velocity at d1
u1_ind = pi*(7000/720)*d1;

% inducer design relative inlet flow velocity
v0_ind = sqrt(cm0_ind.^2 + u0_ind.^2);

sin_beta0 = cm0_ind / v0_ind;
beta0_ind = asin(sin_beta0);
beta0_deg_ind = rad2deg(beta0_ind);

% tangential component of the inducer absolute outlet velocity
cu1_ind = (H_inducer_ft * g) / u1_ind;

% inducer design absolute outflow flow velocity
c1_ind = sqrt(cu1_ind.^2 + cm1_ind.^2);

% inducer design absolute outflow flow angle
alpha1_ind = atan(cm1_ind/cu1_ind);
alpha1_deg_ind = rad2deg(alpha1_ind);

% inducer design relative outlet flow velocity 
v1_ind = sqrt((u1_ind - cu1_ind).^2 + cm1_ind.^2);

% inducer design relative outlet flow angle
beta1_ind = atan(cm1_ind / (u1_ind - cu1_ind));
beta1_deg_ind = rad2deg(beta1_ind);

%%
figure
grid on;
grid minor;
hold on;
axis equal;

pos1 = [0, 0, u0_ind, 0];
quiver2(pos1, 'red');
pos2 = [u0_ind, 0, -v0_ind*cos(beta0_ind), v0_ind*sin(beta0_ind)];
quiver2(pos2, 'blue');
pos3 = [0, 0, cu0_ind, cm0_ind];
quiver2(pos3, 'yellow');

pos4 = [u0_ind, 0, u1_ind, 0];
quiver2(pos4, 'red');
pos5 = [u0_ind+u1_ind, 0, -v1_ind*cos(beta1_ind), v1_ind*sin(beta1_ind)];
quiver2(pos5, 'blue');
pos6 = [u0_ind, 0, cu1_ind, cm1_ind];
quiver2(pos6, 'yellow');

xline(u0_ind, '--');
str1 = {'INLET'};
text(150, -50, str1);
str2 = {'OUTLET'};
text(520, -50, str2);

figure
grid on;
grid minor;
hold on;
q1 = quiver(0, 0, u0_ind, 0, 1);
q1.MaxHeadSize = 0.08;
q2 = quiver(u0_ind, 0, -v0_ind*cos(beta0_ind), v0_ind*sin(beta0_ind), 1);
q2.MaxHeadSize = 0.1;
q3 = quiver(0, 0, cu0_ind, cm0_ind);
q3.MaxHeadSize = 0.1;

q4 = quiver(u0_ind, 0, u1_ind, 0, 1);
q4.MaxHeadSize = 0.08;
q5 = quiver(u0_ind+u1_ind, 0, -v1_ind*cos(beta1_ind), v1_ind*sin(beta1_ind), 1);
q5.MaxHeadSize = 0.1;
q6 = quiver(u0_ind, 0, cu1_ind, cm1_ind, 1);
q6.MaxHeadSize = 0.1;

xline(u0_ind, '--');
str1 = {'INLET'};
text(150, -10, str1);
str2 = {'OUTLET'};
text(525, -10, str2);

legend('u0', 'v0', 'c0', 'u1', 'v1', 'c1', '');
title('Triangolo di velocità INDUCER POMPA LOX');

%%
% inducer inlet tip speed
u0_t_ind = (pi*N)/720 * d_0t;

beta0_t = atan(cm0_ind / u0_t_ind);
beta0_t_deg = rad2deg(beta0_t);

% inducer inlet tip % ipotesi
theta0_t = 9;

% angle of attack at the inlet tip
aoa = theta0_t - beta0_t_deg;

% vane angle at the inducer inlet mean effective diameter d_0
theta0 = atan(d_0t / d0 * tan(deg2rad(theta0_t)));
theta0_deg = rad2deg(theta0);

% vane angle at the inducer inlet hub diameter d_0h
theta0_h = atan(d_0t / d_0h * tan(deg2rad(theta0_t)));
theta0_h_deg = rad2deg(theta0_h);

% inducer inlet flow coefficient
phi_ind = cm0_ind / u0_t_ind;

% theoretical inducer suction specific speed
NSS_ind = 