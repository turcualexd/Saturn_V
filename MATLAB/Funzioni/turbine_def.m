clear; close all; clc;

% dati iniziali (interpolazione dati AIAA + manuali F-1)
c_p       = 2742.2380 ;      % c_p gas generator from interpolation
T_gg      = 1062;  % 
T_turb_e  = 888.38;
p_in      = 63.43; % bar
eps       = 16.4;
p_e       = p_in/eps;
o_f       = 0.416; %
g         = 1.128179;
R         = 294.4417;
mdot      = 77.92;
omega     = 574.7;
U_over_c0 = 0.225;

% ipotesi nostre sy rendimenti (AIAA)
eta_ovrll = 0.6;    % rendimento totale turbina (termodinamic + meccanico)
noz_ar    = 9.7;    % aspect ratio nozzle
k_n       = 0.96;   % perdita spouting velocity
eta_n     = k_n^2;  % rendimento ugello
eps_nt    = 0.97;   %area coefficient throat
eps_ne    = 0.95;   %area cofficient exit
k_blade   = 0.89;   % rotor and stator blade velocity coefficient
eps_blade = 0.95;   % stator and rotor blade exit area coefficient


% struttura turbina 

% nozzle inlet - nozzle out --> rotor --> stator --> rotor

%% Inlet nozzle -> suppongo che il 94% del salto di dh avvenga nei nozzle

dh_tot     = c_p*T_gg * (1 - (1/eps)^((g-1)/g)); % salto totale turbina.
dh_nz      = (1 - 0.06)*dh_tot;
p_e_nzz    = p_in*(1 - dh_nz/(c_p*T_gg))^(g/(g-1));
C_zero     = sqrt(2*dh_nz);
C_1        = k_n * C_zero; 
U          = C_zero*U_over_c0;

q_n_re      = ((1 - k_n^2)*C_1^2)/(k_n^2);
T_e_nzz_id  = T_gg - dh_nz/c_p;
T_e_nzz_re  = T_e_nzz_id + q_n_re/c_p;

%relazione per efficienza blade massima
alfa_1     = acos(4*U/C_1);
alfa_1_deg = alfa_1*(180/pi);

beta_1     = atan(C_1*sin(alfa_1)/(C_1*cos(alfa_1) - U));
beta_1_deg = rad2deg(beta_1);

V_1        = C_1*sin(alfa_1)/sin(beta_1);


%% solving system 
% ipotesi --> mach uscita 0.4, alfa2 = alfa3, c3 = kc2, v2= kv1 e v4 = kv3
% risolvendo il sistema troviamo alfa 2 e C2

x0     = [C_1, pi/3];
x_sol  = fsolve(sist, x0); % troviamo rispettivamente C_2 e alfa2
C_2    = x_sol(1);
alfa_2 = x_sol(2);

%% first rotor exit - stator  inlet

U_vec       = U*[1, 0];
T_e_rot1_id = T_e_nzz - dh_rot1/c_p;
q_rot1_re   = (1 - k_blade^2)*(V_1^2 / 2) + (1 - k_n^2)*dh_rot1;
T_e_rot1_re = T_e_rot1_id + q_nr/c_p;
C_2_vec     = [C_2*cos(pi + alfa_2), C_2*sin(pi + alfa_2)];
V_2_vec     = C_2_vec - U_vec;
p_e_rot1    = p_e_nzz*(1 - dh_nz/(c_p*T_e_nzz_re))^(g/(g-1));
beta_2      = acos(V_2_vec(1)/norm(V_2_vec));

%% stator outlet - 2nd rotor inlet

alfa_3      = alfa_2;
C_3         = k_blade*C_2;
C_3_vec     = C_3*[cos(-alfa_3), sin(-alfa_3)];
V_3_vec     = C_3_vec - U_vec;
q_bs        = (1 - k_blade^2)*(C_2^2 / 2) + (1 - k_n^2)*dh_stat;
T_e_stat_id = T_e_rot1_re - dh_stat/c_p;
T_e_stat_re = T_e_stat_id + q_bs/c_p;
p_e_stat    = p_e_rot1*(1 - dh_stat/(c_p*T_e_rot1_re))^(g/(g-1));


%% 2nd rotor outlet

alfa_4      = pi/2;
C_4         = 0.4*sqrt(g*R*T_turb_e);
C_4_vec     = C_4 * [cos(-alfa_4), sin(-alfa_4)];
V_4_vec     = C_4_vec - U_vec;
q_rot2_re   = (1 - k_blade^2)*(norm(V_3_vec)^2 / 2) + (1 - eta_n)*dh_rot1;
T_e_rot2_id = T_e_stat_re - dh_rot1/c_p;
T_e_rot2_re = T_e_rot2_id + q_rot2_re/c_p;

%% diagrammi

figure;
q1 = quiver(0,0,C_1*cos(alfa_1), C_1*sin(-alfa_1),1,'r');
q1.MaxHeadSize = 0.05;
hold on;
q2 = quiver(0,0, V_1*cos(-beta_1), V_1*sin(-beta_1),1,'r');
q2.MaxHeadSize = 0.09;
q3 = quiver(V_1*cos(-beta_1), V_1*sin(-beta_1), 1e-3*U, 0, 1, 'b');
q2.MaxHeadSize = 0.08;
grid on
grid minor
plot(linspace(-0.5, 1.5, 100), C_1_vec(2)*ones(100,1), ':r', 'LineWidth',1);

q4 = quiver(V_1_vec(1), V_1_vec(2), C_2_vec(1), C_2_vec(2), 1, 'r');
q4.MaxHeadSize = 0.09;
q5 = quiver(V_1_vec(1), V_1_vec(2), V_2_vec(1), V_2_vec(2), 1, 'r');
q5.MaxHeadSize = 0.09;
q6 = quiver(V_1_vec(1) + V_2_vec(1), V_1_vec(2) + V_2_vec(2), U, 0, 1, 'b');
q6.MaxHeadSize = 0.08;

q7 = quiver(V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), C_3_vec(1), C_3_vec(2), 1, 'r');
q7.MaxHeadSize = 0.09;
q8 = quiver(V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), V_3_vec(1), V_3_vec(2), 1, 'r');
q8.MaxHeadSize = 0.09;
q9 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), U, 0, 1, 'b');
q9.MaxHeadSize = 0.08;

q10 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), C_4_vec(1), C_4_vec(2), 1, 'r');
q10.MaxHeadSize = 0.09;
q11 =quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), V_4_vec(1), V_4_vec(2), 1, 'r');
q11.MaxHeadSize = 0.09;
q12 = quiver(V_1_vec(1) + C_2_vec(1) + V_3_vec(1) + V_4_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2) + V_4_vec(2), U, 0, 1, 'b');
q12.MaxHeadSize = 0.08;

plot(linspace(-0.5, 1.5, 100),(C_1_vec(2) + C_2_vec(2))*ones(100,1), ':r', 'LineWidth',1);
plot(linspace(-0.5, 1.5, 100), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2))*ones(100,1), ':r', 'LineWidth',1);
plot(linspace(-0.5, 1.5, 100), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2) + C_4_vec(2))*ones(100,1), ':r', 'LineWidth',1);



