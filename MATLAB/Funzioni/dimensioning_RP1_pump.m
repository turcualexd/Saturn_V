clear; close all; clc;

%RP1 pump

%design parameters
eta_p  = 0.760;
P01    = 45;         %psia
P02    = 1860;       %psia
deltaP = P02 - P01;  %psia
dH     = 5175.445;   % head rise ft
Rtip   = 0.59436/2;  %m raggio impeller punta
Dtip   = Rtip*2;
omega  = 575.12;     %rad/s
Q      = 0.9617;      %m3s
N      = 5488;       %rpm
g      = 9.81;       % 
theta  = deg2rad(10);% rad: angolo di leading edge inducer
rD     = 0.81;       %rapporto diametro hub/tip
rho_RP1= 808.9324;   %kg/m3
rho_USA= 50.5;       %lbm/ft3
beta2  = deg2rad(25);% rad: back angle impeller at outlet
b2     = 0.04318;    % m spessore paletta alla punta
m_dot  = 796.5082;
P1     = 310264;     %Pa 
Psat   = 137895;     %Pa

%psi  pump head coefficient
%NPSP (uguale a NPSH)
%sigma parametro thoma
%Cm velocitò linea meridionale
NPSH = 55; %feet (critical)
NPSP = NPSH/(32.2 * rho_USA); %psia

[dH_m, ~ ,psi, pwr, ds, Ns, prod,~, ~, ~, ~, ~, Nss, ~, ~] = ...
turbopump(eta_p, P01, P02, rho_USA, Q, m_dot, Dtip, b2, beta2, ...
    omega, P1, Psat, rho_RP1, rD, theta);
Q_gpm = Q/(6.309*1e-5);
%hp
Nss_imp  = 11000; 
NPSH_imp = ((N*(Q_gpm)^0.5)/(Nss_imp))^1.33;
dH_ind   = NPSH_imp - NPSH; 

%hp psi inducer
psi_ind = 0.06;
U_tip_imp = omega*Rtip;
U_tip_ind = sqrt(dH_ind*32.2/psi_ind)/3.281; %m/s

d_tip_ind = (720/(pi*N))*(U_tip_ind*3.281); % inches
l_ind     = 0.4 * d_tip_ind;                 % inches


