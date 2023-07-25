clear
close all
clc

%% DATI POMPA
% Inducer
P1 = 310264;          % Pa
Psat = 137895;        % Pa
pho = 808.9324;       % kg/m^3
rapportoD = 0.81;

% Impeller
eta_p = 0.760;
P01 = 45;             % psia
P02 = 1860;           % psia
pho1 = 50.5;          % lbm/ft^3
Q = 0.9617;           % m^3/s
m_dot = 796.5082;     % kg/s
Dt2 = 0.59436;        % m
b2 = 0.04318;         % m
beta2 = deg2rad(25);
omega = 574.020;      % rad/s

U2 = omega*(Dt2/2);
Dt1 = rapportoD * Dt2;
R1 = Dt1/2;
U1 = omega*R1;

[H_m_RP1, delta_h_RP1, psi_RP1, P_RP1, ds_RP1, ns_RP1, prod_RP1] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, omega);

% Le ipotesi sottostanti sono basate su considerazioni tratte da:
% author = {D.K. Huzel, D.H. Huang},
% title  = {Design of Liquid Propellant Rocket Engines, vol. 2},
% year   = {1967},
% publisher = {Scientifical and Technichal Division NASA}

%% 1° tratto: INDUCER
% gravitational constant
g = 32.2;           % ft/s
N = 5488;           % rpm
Q = 15250;          % gpm
NPSH_c_RP1 = 55;    % ft
H = (144.*(P02 - P01) ./ pho1); % ft

NS = (N*sqrt(Q)) ./ (H^0.75);

% Parametro di Thoma
thoma = NPSH_c_RP1 / H;
NSS = NS / (thoma)^0.75;

% Ipotesi pg. 224
NSS_impeller = 11000; 

% NPSH impeller
NPSH_imp = (N*sqrt(Q) ./ NSS_impeller).^(4/3);

H_inducer_ft = NPSH_imp - NPSH_c_RP1;
perc_ind = (H_inducer_ft / H) * 100;

% Ipotesi pg. 224
psi_ind = 0.06; 

u_tip_mean = sqrt((H_inducer_ft*32.2)/psi_ind);

% Inducer mean tip diameter
d_t = (720/(pi*N)) .* u_tip_mean; % in

L_ind = d_t * 0.4;
% tip diameter at the inducer inlet
d_0t = d_t + 2*L_ind/2 * tan(deg2rad(7)); % in
% tip diameter at the inducer outlet
d_1t = d_t - 2*L_ind/2 * tan(deg2rad(7)); % in

% Ipotesi pg. 224
rd = 0.3;

% Mean hub diameter
d_h = d_t * rd; % in
% Hub diameter at the inducer inlet
d_0h = d_h - 2*L_ind/2 * tan(deg2rad(14)); % in
% Hub diameter at the inducer outlet
d_1h = d_h + 2*L_ind/2 * tan(deg2rad(14)); % in

Q_inducer = Q * (1 + 0.032 + 0.0175); % gpm (pg. 225)

cm0_ind = Q_inducer/ (3.12 * (pi/4) * (d_0t.^2 - d_0h.^2)); % in ft/s
cm0_ind_SI = cm0_ind / 3.281; % in m/s
% Assunzione pg. 226
cu0_ind = 0; % in m/s
cm1_ind = Q_inducer/ (3.12 * (pi/4) * (d_1t.^2 - d_1h.^2)); % in ft/s
cm1_ind_SI = cm1_ind / 3.281; % in m/s

% Inducer mean effective diameter inlet
d0 = sqrt((d_0t.^2 + d_0h.^2) / 2); % in

% Inducer peripheral velocity at d0
u0_ind = pi*(N/720)*d0; % in ft/sec
u0_ind_SI = u0_ind / 3.281; % in m/s

% Inducer mean effective diameter outlet
d1 = sqrt((d_1t.^2 + d_1h.^2) / 2); % in

% Inducer peripheral velocity at d1
u1_ind = pi*(N/720)*d1; % in ft/sec
u1_ind_SI = u1_ind / 3.281; % in m/s

% Inducer design relative inlet flow velocity
v0_ind = sqrt(cm0_ind.^2 + u0_ind.^2); % in ft/s
v0_ind_SI = v0_ind / 3.281; % in m/s

sin_beta0 = cm0_ind / v0_ind;
beta0_ind = asin(sin_beta0);
beta0_deg_ind = rad2deg(beta0_ind);

% Tangential component of the inducer absolute outlet velocity
cu1_ind = (H_inducer_ft * g) / u1_ind; % in ft/s
cu1_ind_SI = cu1_ind / 3.281; % in m/s

% Inducer design absolute outflow flow velocity
c1_ind = sqrt(cu1_ind.^2 + cm1_ind.^2); % in ft/s
c1_ind_SI = c1_ind / 3.281; % in m/s

% Inducer design absolute outflow flow angle
alpha1_ind = atan(cm1_ind/cu1_ind);
alpha1_deg_ind = rad2deg(alpha1_ind);

% Inducer design relative outlet flow velocity 
v1_ind = sqrt((u1_ind - cu1_ind).^2 + cm1_ind.^2); % in ft/s
v1_ind_SI = v1_ind / 3.281; % in m/s

% Inducer design relative outlet flow angle
beta1_ind = atan(cm1_ind / (u1_ind - cu1_ind));
beta1_deg_ind = rad2deg(beta1_ind);

% Velocity diagram INDUCER: plot
figure
grid on;
grid minor;
hold on;
axis equal;

pos1_ind = [0, 0, u0_ind_SI, 0];
quiver2(pos1_ind, 'blue');
pos2_ind = [u0_ind_SI, 0, -v0_ind_SI*cos(beta0_ind), v0_ind_SI*sin(beta0_ind)];
quiver2(pos2_ind, 'blue');
pos3_ind = [0, 0, cu0_ind, cm0_ind_SI];
quiver2(pos3_ind, '#00FFFF');

pos4_ind = [u0_ind_SI, 0, u1_ind_SI, 0];
quiver2(pos4_ind, 'blue');
pos5_ind = [u0_ind_SI+u1_ind_SI, 0, -v1_ind_SI*cos(beta1_ind), v1_ind_SI*sin(beta1_ind)];
quiver2(pos5_ind, 'blue');
pos6_ind = [u0_ind_SI, 0, cu1_ind_SI, cm1_ind_SI];
quiver2(pos6_ind, '#00FFFF');

pos7_ind = [u0_ind_SI+u1_ind_SI-v1_ind_SI*cos(beta1_ind), 0, 0, cm1_ind_SI];
quiver2(pos7_ind, 'k');

xline(u0_ind_SI, '--');
str1 = {'INLET'};
text(25, -20, str1);
str2 = {'OUTLET'};
text(105, -20, str2);
title('Diagramma di velocità INDUCER POMPA RP-1')
xlabel('velocità tangenziale [m/s]')
ylabel('velocità meridionale [m/s]')

% Inducer inlet tip speed
u0_t_ind = (pi*N)/720 * d_0t; % in ft/sec
u0_t_ind_SI = u0_t_ind / 3.281; % in m/s

beta0_t = atan(cm0_ind / u0_t_ind);
beta0_t_deg = rad2deg(beta0_t);

% Inducer inlet tip % ipotesi pg. 226
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
NSS_ind = ((8150*(1-2*phi_ind.^2)).^(0.75))/(phi_ind) .* (1-rd.^2).^0.5;

% vane angle 14.30° ipotesi pg. 226
theta1 = 14.30;

% vane angle at the inducer outlet tip diameter d_1t
theta1_t = atan(d1 / d_1t * tan(deg2rad(theta1)));
theta1_t_deg = rad2deg(theta1_t);

% vane angle at the inducer outlet hub diameter d_1h
theta1_h = atan(d1 / d_1h * tan(deg2rad(theta1)));
theta1_h_deg = rad2deg(theta1_h);

% z = 3 vanes pg. 226
z = 3;
% vane pitch at the mean tip diameter dt
P_i = (pi*d_t) / z; % in

% Chord length at vane tip
C_i = L_ind / sin((deg2rad(theta0_t) + theta1_t)/2); % in

% Inducer solidity
S_v = C_i / P_i;

%% 2° tratto: IMPELLER
% Ipotesi pg. 216
H_e = 0.30*H;
H_imp = H - H_e;

Dt2_f = Dt2/0.3048; % in ft
Dt2_inches = Dt2 * 39.3701;
R2 = Dt2_f/2;
Dt1_f = rapportoD * Dt2_f; % in ft
R1 = Dt1_f/2;

% Area normal to the radial flow at the impeller inlet
A1 = pi*(R1.^2); % ft^2

% Area normal to the radial flow at the impeller outlet
A2 = pi*(R2.^2-R1.^2); % ft^2

% Radial/Meridional component of the absolute inlet flow velocity
cm1 = Q / (448.8*A1); % in ft/s
cm1_SI = cm1 / 3.281; % in m/s

% Radial/Meridional component of the absolute inlet flow velocity
% ipotesi da Sistemi Energetici 
% DESIGN OF A CENTRIFUGAL PUMP FOR AN EXPANDER CYCLE ROCKET ENGINE, 
% Prof. Ing. Aldo BONFIGLIOLI
cm2 = 1.5*cm1; % in ft/s
cm2_SI = cm2 / 3.281; % in m/s

% Impeller peripheral velocity at inlet
u1 = U1 * 3.281; % in ft/s
u1_SI = u1 / 3.281; % in m/s

% Impeller peripheral velocity at outlet
u2 = U2 * 3.281; % in ft/s
u2_SI = u2 / 3.281; % in m/s

% Tangential component of the absolute outlet flow velocity
cu2 = u2 - (cm2 / tan(beta2)); % in ft/s
cu2_SI = cu2 / 3.281; % in m/s

% Tangential component of the absolute inlet flow velocity
cu1 = (u2*cu2 - H*g) / u1; % in m/s;
cu1_SI = cu1 / 3.281; % in m/s;

% Absolute inlet velocity of the flow
c1 = sqrt(cu1.^2 + cm1.^2); % in ft/s
c1_SI = c1 / 3.281; % in m/s

% Absolute outlet velocity of the flow
c2 = sqrt(cu2.^2 + cm2.^2); % in ft/s
c2_SI = c2 / 3.281; % in m/s

vu1 = sqrt((u1 - cu1).^2 + cm1.^2); % in ft/s
vu1_SI = vu1 / 3.281; % in m/s
vm1 = cm1;
vm1_SI = vm1 / 3.281; % in m/s

% Inlet flow velocity relative to the impeller
v1 = sqrt(vu1.^2 + vm1.^2); % in ft/s
v1_SI = v1 / 3.281; % in m/s

vu2 = sqrt((u2 - cu2).^2 + cm2.^2); % in ft/s
vu2_SI = vu2 / 3.281; % in m/s
vm2 = cm2;
vm2_SI = vm2 / 3.281; % in m/s

% Outlet flow velocity relative to the impeller
v2 = sqrt(vu2.^2 + vm2.^2); % in ft/s
v2_SI = v2 / 3.281; % in m/s

alpha1 = atan(cm1/cu1);
alpha1_deg = rad2deg(alpha1);

% Impeller inlet vane angle
beta1 = atan(vm1/vu1);
beta1_deg = rad2deg(beta1);

alpha2 = atan(cm2/cu2);
alpha2_deg = rad2deg(alpha2);

% Impeller discharge vane angle
% beta2_controllo = atan(vm2/vu2);
% beta2_deg_controllo = rad2deg(beta2_controllo);

% Velocity diagram impeller: plot
figure
grid on;
grid minor;
hold on;
axis equal;

pos1 = [0, 0, u1_SI, 0];
quiver2(pos1, 'blue');
pos2 = [u1_SI, 0, -vu1_SI, vm1_SI];
quiver2(pos2, 'blue');
pos3 = [0, 0, u1_SI-vu1_SI, cm1_SI];
quiver2(pos3, '#00FFFF');

pos4 = [u1_SI, 0, u2_SI, 0];
quiver2(pos4, 'blue');
pos5 = [u1_SI+u2_SI, 0, -vu2_SI, vm2_SI];
quiver2(pos5, 'blue');
pos6 = [u1_SI, 0, u2_SI-vu2_SI, cm2_SI];
quiver2(pos6, '#00FFFF');

pos7 = [u1_SI-vu1_SI, 0, 0, cm1_SI];
quiver2(pos7, 'k');
pos8 = [u1_SI+u2_SI-vu2_SI, 0, 0, cm2_SI];
quiver2(pos8, 'k');

xline(u1_SI, '--');
str1 = {'INLET'};
text(60, -20, str1);
str2 = {'OUTLET'};
text(210, -20, str2);
title('Diagramma di velocità IMPELLER POMPA RP-1')
xlabel('velocità tangenziale [m/s]')
ylabel('velocità meridionale [m/s]')

%% 3° tratto: Pump casing
% Design of a double-volute (spaced 180°) single-discharge-type casing
% double-volute: si divide il flusso in 2 correnti uguali usando due 
% linguette distanziate di 180°

% Ipotesi pg. 231
Kv = 0.337; % design factor

% Average flow velocity in the volute
c3 = Kv * sqrt(2*g*H); % in ft/s
c3_SI = c3 / 3.281; % in m/s

% area of a volute section at a given angular location \theta_volute from the tongue
a_thetav_const = Q / (3.12 * 360 * c3);

a_volute_45 = a_thetav_const * 45; % in^2
a_volute_90 = a_thetav_const * 90; % in^2
a_volute_135 = a_thetav_const * 135; % in^2
a_volute_180 = a_thetav_const * 180; % in^2

thetav = [45:0.01:180];
a_volute_theta = [];
for i = 1:length(thetav)
    a_v = a_thetav_const*thetav(i);
    a_volute_theta = [a_volute_theta; a_v];
end

% area of the volute throat section
a_v = 2*a_volute_180;

% volute angle is similar to alpha_2_deg = 14.7757° 
% for geometrical construction pg. 231
alpha_v_deg = 15;
alpha_v = deg2rad(alpha_v_deg);

% radius at which start the volute tongue
% ipotesi: si assume un gioco del 5% pg. 231
r_t = (Dt2_inches / 2)*1.05; % in

% width at the bottom of the trapezoidal volute section
b3_SI = 1.75*b2; % in m
b3 = b3_SI *39.37; % in

% pitch diameter of the diffuser throats
d3 = Dt2_inches*(c2/c3); % in

% Ipotesi pg. 231
% transizione dalla forma del diffusore ad un cerchio: d_transition
% taper angle di 10° e lunghezza nozzle di 10 inches

% da costruzione grafica
a_transition = a_v + b3; %in^2 !!!
d_transition = sqrt((4*a_transition) / pi);

% diameter of the discharge nozzle
d_e = d_transition + 2*10*tan(deg2rad(5));
A_e = pi*(d_e^2)/4;

% flow velocity at the nozzle inlet
% v_nozzle_i = Q / (3.12*30.68); % in ft/s
v_nozzle_i = c2; % assunzione
v_nozzle_i_SI = v_nozzle_i / 3.281; % in m/s

% flow velocity at the nozzle exit
v_nozzle_e = Q / (3.12*A_e); % in ft/s
v_nozzle_e_SI = v_nozzle_e / 3.281; % in m/s
