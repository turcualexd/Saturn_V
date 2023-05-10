clear; close all; clc;

% dati iniziali (interpolazione dati AIAA + manuali F-1)
c_p       = 2742.2380 ;      % c_p gas generator from interpolation
T_gg      = 1062;  % 
p_in      = 63.43; % bar
eps       = 16.4;
p_e       = p_in/eps;
o_f       = 0.416; %
g         = 1.128179;
mdot      = 77.92;
omega     = 574.7;
U_over_c0 = 0.2;

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
T_e_nzz     = T_e_nzz_id + q_n_re/c_p;

%relazione per efficienza blade massima
alfa_1     = acos(4*U/C_1);
alfa_1_deg = alfa_1*(180/pi);

beta_1     = atan(C_1*sin(alfa_1)/(C_1*cos(alfa_1) - U));
beta_1_deg = rad2deg(beta_1);

V_1        = C_1*sin(alfa_1)/sin(beta_1);

%% first rotor exit - stator  inlet

dh_rot1     = dh_tot*0.06/3;
V_2         = sqrt(k_blade^2 * V_1^2 + 2*dh_rot1);
q_blade_rh  = (1 - k_blade^2)*(V_1^2 / 2) + (1 - k_n^2)*dh_rot1;
T_e_rot1_id = T_e_nzz - dh_rot1/c_p;
T_e_rot1    = T_e_rot1_id + q_blade_rh/c_p;

%% statore

dh_stat     = dh_rot1;




%% procedimento per ricavare gli altri angoli 
% ipotizziamo che l'angolo asosluto in sucita sia di alfa4 = 90°, inoltre
% l'angolo di uscita dallo statore è uguale all'angolo di ingresso. 

%fun = @solve_sistem;
x0  = [pi/3 C_1];
f   = @solve_sistem;
x   = fsolve(f, x0);




