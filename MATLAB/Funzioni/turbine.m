clear; close all; clc;

%interpolazione dati
%..

c_p     = 2742.2380 ;      % c_p gas generator from interpolation
T_gg    = 1062;  % K
p_in    = 63.43; % bar
eps     = 16.4;
p_e     = p_in/eps;
o_f     = 0.416; %
g       = 1.128179;
C_zero  = sqrt(2*c_p*T_gg * (1 - (1/eps)^((g-1)/g) ));
mdot    = 77.92;
omega   = 574.7;
eta     = 0.605;
pwr     = 39.45*1e6; %watt

%pitch line velocity (omega * R_turbina)
U         = 256; %m/s
r_turb    = U/omega; %raggio turbina-rotore (torna con manuale nasa)
vel_ratio = U/C_zero;
psi       = 0.5*(C_zero/U)^2;

% ipotesi nostre - nozzle turbina
k_n = 0.89; % valore più basso dato da AIAA per calcolare c_1 (velocità 
            % assoluta in uscita dai nozzle)
eta_n  = k_n^2;
 eta_m = eta/eta_n;
C_1 = k_n * C_zero;
q_n = ((1 - k_n^2) * C_1^2) / (k_n); % nozzle reheat
eps_n = 0.95; %from AIAA

alfa_1     = acos(4*U/C_1);
alfa_1_deg = alfa_1*(180/pi);

% diagrammi di velocità (km/s)
% criteri di progetto: vogliamo che il vettore uscente dalla turbina sia
% piccolo, e assiale (90 gradi), assumiamo che ogni stadio abbia una
% perdita dovuta a attrito di k = ...

figure;
q1 = quiver(0,0,1e-3*C_1*cos(alfa_1), 1e-3*C_1*sin(-alfa_1),1,'r');
q1.MaxHeadSize = 0.05;
%gtext('U_1')
hold on
C_1_vec =1e-3 * [C_1*cos(-alfa_1), C_1*sin(-alfa_1)];
U_vec   = 1e-3*[U, 0];
V_1_vec = C_1_vec - U_vec;
V_1     = norm(V_1_vec);
beta_1  = acos(V_1_vec(1,1)/V_1);
beta_1_deg = beta_1*180/pi;
q2 = quiver(0,0, V_1*cos(-beta_1), V_1*sin(-beta_1),1,'r');
q2.MaxHeadSize = 0.09;
q3 = quiver(V_1*cos(-beta_1), V_1*sin(-beta_1), 1e-3*U, 0, 1, 'b');
q2.MaxHeadSize = 0.08;
grid on
grid minor
plot(linspace(-0.5, 1.5, 100), C_1_vec(2)*ones(100,1), ':r', 'LineWidth',1);

eta_vel = 0.7;

% V_2 = eta_vel * V_1;
% V_3 = eta_vel * V_2;
% V_4 = eta_vel * V_3;

V_2 = eta_vel * V_1;
%C_3 = eta_vel * C_2;
%V_4 = eta_vel * V_3;

%vincolo uscita
alfa_4 = 79*pi/180;

dh_iso_tot = c_p * T_gg * (1 - (1/eps)^((g-1)/g));
dh_real    = pwr/mdot;

%e_tot = U*(C_1*cos(alfa_1) + C_2*cos(alfa_2) + C_3*cos(alfa_3) +
%C_4*cos(alfa_4)) con alfa_3 = alfa_1 , alfa_2 incognita e alfa_4 = 90;
alfa_3 = alfa_1;
beta_3 = beta_1;

%alfa_2 = acos(dh_real/(U*C_2) - (C_1*cos(alfa_1) + C_3*cos(alfa_3) + C_4*cos(alfa_4))/C_2 );
beta_2 = acos((1/(eta_vel*1e3*V_1)) *  (dh_real/U - V_1*1e3*cos(beta_1) - 0.9*V_2*1e3*cos(beta_3) - U));

% C_2_vec = 1e-3*[C_2*cos(alfa_2 + pi), C_2*sin(alfa_2 + pi)];
% V_2_vec = C_2_vec - U_vec;
V_2_vec = [V_2*cos(beta_2 + pi), V_2*sin(beta_2 + pi)];
C_2_vec = V_2_vec + U_vec;

% C_3_vec = 1e-3*[C_3*cos(-alfa_3), C_3*sin(-alfa_3)];
% V_3_vec = C_3_vec - U_vec;
V_3 = 0.9*V_2;
V_3_vec = [V_3*cos(-beta_3), V_3*sin(-beta_3)];
C_3_vec = V_3_vec + U_vec;

V_4 = eta_vel * norm(V_3_vec);
beta_4 = acos(1e-3*U/(V_4));
% C_4_vec = 1e-3*[C_4*cos(pi + alfa_4), C_4*sin(pi + alfa_4)];
% V_4_vec = C_4_vec - U_vec;
V_4_vec = [V_4*cos(pi + beta_4), V_4*sin(pi + beta_4)];
C_4_vec = V_4_vec + U_vec;


q4 = quiver(V_1_vec(1), V_1_vec(2), C_2_vec(1), C_2_vec(2), 1, 'r');
q4.MaxHeadSize = 0.09;
q5 = quiver(V_1_vec(1), V_1_vec(2), V_2_vec(1), V_2_vec(2), 1, 'r');
q5.MaxHeadSize = 0.09;
q6 = quiver(V_1_vec(1) + V_2_vec(1), V_1_vec(2) + V_2_vec(2), 1e-3*U, 0, 1, 'b');
q6.MaxHeadSize = 0.08;

q7 = quiver(V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), C_3_vec(1), C_3_vec(2), 1, 'r');
q7.MaxHeadSize = 0.09;
q8 = quiver(V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), V_3_vec(1), V_3_vec(2), 1, 'r');
q8.MaxHeadSize = 0.09;
q9 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), 1e-3*U, 0, 1, 'b');
q9.MaxHeadSize = 0.08;

q10 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), C_4_vec(1), C_4_vec(2), 1, 'r');
q10.MaxHeadSize = 0.09;
q11 =quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), V_4_vec(1), V_4_vec(2), 1, 'r');
q11.MaxHeadSize = 0.09;
q12 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1) + V_4_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2) + V_4_vec(2), 1e-3*U, 0, 1, 'b');
q12.MaxHeadSize = 0.08;

plot(linspace(-0.5, 1.5, 100),(C_1_vec(2) + C_2_vec(2))*ones(100,1), ':r', 'LineWidth',1);
plot(linspace(-0.5, 1.5, 100), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2))*ones(100,1), ':r', 'LineWidth',1);
plot(linspace(-0.5, 1.5, 100), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2) + C_4_vec(2))*ones(100,1), ':r', 'LineWidth',1);


