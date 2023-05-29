%% CALCOLO DELLO SCAMBIO TERMICO NELL'UGELLO AGGIUNTIVO
% scambio convettivo, liquid rocket design da pagina 110

clear
clc

% valori nel punto di aggiunta dell'estensione dell'ugello

T_c = 6429.6;     % R            temperatura in camera di combustione
p_c = 1124.85;    % psia         

g   = 1.2439;     % /            Gamma dei gas combusti 
mm  = 0.4983;     % lb/mol       Massa molare 
cp  = 0.448;      % btu/lb*F     Specific heat at constant pressure a 1:10
c_star = 5929.2;  % ft/s         

T_g = 3012.64;    % R            Temperatura gas combusti in 1:10
R_t = 17.5;       % in           Raggio in gola

f_aw   = 0.94;    % /            Turbolent boundary layer recovery factor [val medio 0.9-0.98]
f_wg   = 0.8;     % T_wg/Tc 
R_d = 2001.8;       % in^2*s*F/btu Thermal resistance caused by the solid deposition

R = 6.96;         % in           Radius of curvature of nozzle contour (parte nozzle)

rapp_A  = 1/10;    % /           Rapporto area gola / inizio ugello aggiunto

sigma = 0.7;      %              Correction factor for property variations across the boundary layer (all'uscita dal nozzle)

gg = 32.2;        % ft/s^2       Accelerazione gravitazionale


% OSS: R_d, f_aw, sigma sono presi da pagina 111 del liquid rocket design

% CALCOLI:

T_wg = f_wg * T_c;     % R     Hot-gas-side local chamber-wall temperature
T_aw = f_aw * T_c;      % R     Adiabatic wall temperature of the gas


mu = 46.6e-10 * mm^0.5 * T_g^0.6;       % viscosità
Pr = 4*g / ( 9*g -5 );                  % considero costante 

hg = (0.026/((R_t*2)^0.2)) * (mu^0.2 * cp / (Pr^0.6)) * (p_c * gg / c_star)^0.8 * (((R_t*2) / R)^0.1 ) * (rapp_A)^0.9 * sigma; % gas side heat transfer coefficient [btu/in^2*sec*F]
h_gc = 1 / ( (1/hg) + R_d ); % btu/in^2*sec*F     overall gas-side thermal conductance 

q = h_gc * (T_aw - T_wg); %Heat flux or heat transferred across the stagnant gas film per unit surface area per unit time [btu/in^2*sec]

q_real = q * 1055.06/ 0.00064516    % Heat flux [J/m^2*s]

%%


