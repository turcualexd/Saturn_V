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
U      = 256; %m/s
r_turb = U/omega; %raggio turbina-rotore (torna con manuale nasa)
vel_ratio = U/C_zero;
psi = 0.5*(C_zero/U)^2;

% ipotesi nostre - nozzle turbina
k_n = 0.89; % valore più basso dato da AIAA per calcolare c_1 (velocità 
            % assoluta in uscita dai nozzle)
eta_n  =k_n^2;
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
quiver(0,0,1e-3*C_1*cos(alfa_1), 1e-3*C_1*sin(-alfa_1),'r');
hold on
C_1_vec =1e-3 * [C_1*cos(-alfa_1), C_1*sin(-alfa_1)];
U_vec   = 1e-3*[U, 0];
V_1_vec = C_1_vec - U_vec;
V_1     = norm(V_1_vec);
beta_1  = acos(V_1_vec(1,1)/V_1);
beta_1_deg = beta_1*180/pi;
quiver(0,0, V_1*cos(-beta_1), V_1*sin(-beta_1), 'r');
quiver(V_1*cos(-beta_1), V_1*sin(-beta_1), 1e-3*U, 0, 'b');
grid on





dh_iso_tot = c_p * T_gg * (1 - (1/eps)^((g-1)/g));
