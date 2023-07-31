clear; clc;
%% Combustion chamber

L_star = 1; % lunghezza caratteristica, m
R_t = 0.4445; % raggio di gola, m 
A_t = pi*R_t^2; % area di gola, m^2
V_cc = L_star * A_t; % volume camera di combustione, m^3
% volume compreso il convergente (ciclindro + convergente)

% calcolo c_star
R = 8314;  % [kJ / mol * K]
g0 = 9.81;  % [m/s^2]
gamma = 1.1777; 
GAMMA = sqrt(gamma.*(2/(gamma+1)).^((gamma+1)/(gamma-1)));% Vandencherkove function
T_c = 3572 ; % temperatura in camera; K
p_c = 77.6e5; % pressione in camera, Pa
MM = 22.2095; % massa molare, Kg/kmol 

c_star = (p_c .* A_t) ./ (GAMMA .* (p_c ./ (sqrt((R / MM).*T_c))) .* A_t); % velocità caratteristica, m/s
t_r = L_star/c_star; % tempo residenza, s
%% misure della parte cilindrica

eps = 1.307; % A_cil/A_t rateo di contrazione

% V_c è il volume con il convergente, assumiamo 
% il convergente come una % del volume della porzione
% della camera cilindrica

A_cil = eps*A_t;
R_cil = sqrt(A_cil/pi);

% studio parte convergente
a = (R_cil-R_t); 
theta_deg = 13; % deg
theta = deg2rad(theta_deg);

L_conv = a/tan(theta);
V_conv = ((A_cil + A_t + sqrt(A_cil*A_t))*L_conv)/3; %volume parte convergente

% ritorniamo al cilindro
L_cil = (V_cc -V_conv)/(A_cil);

V_cil = L_cil*A_cil; % volume corretto della parte cilindrica

%% area interna delle pareti

A_tot = 2*L_cil*sqrt(pi*eps*A_t) + csc(theta)*(eps-1)*A_t; % area totale della camera di combustione (cilindro + conv)

% apo = sqrt(L_conv^2 + a^2)
% A_conv = pi*(R_c +R_t)*apo;
% A_cil = 2*pi*R_c*L_c;
% A_tot2 = A_conv + A_cil;