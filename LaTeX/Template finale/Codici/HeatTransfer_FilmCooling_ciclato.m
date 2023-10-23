%% CALCOLO DELLO SCAMBIO TERMICO NELL'UGELLO AGGIUNTIVO
% scambio convettivo, liquid rocket design da pagina 110

clear
clc
close all

% valori nel punto di aggiunta dell'estensione dell'ugello

T_c_teorica = 6429.6;       % R            temperatura in camera di combustione
T_c = T_c_teorica*0.975^2;    % R            fattore correttivo di c* pari a 0.975 
p_c = 1124.85;              % psia         

g   = 1.1754;               % /            Gamma dei gas combusti 
mm  = 22.26;                % lb/mol       Massa molare gas combusti (RPA)
cp  = 1.1225;               % btu/lb*F     Specific heat at constant pressure a 1:10
c_star = 5929.2;            % ft/s         

R_t = 17.5;                 % in           Raggio di gola

f_aw   = 0.93;              % /            Turbolent boundary layer recovery factor [val medio 0.9-0.98]
f_wg   = 0.8;               %              T_wg/Tc 
R_d = 2001.8;               % in^2*s*F/btu Thermal resistance caused by the solid deposition

R = (27.38+6.97)/2;         % in           Radius of curvature of nozzle contour (parte nozzle)

rapp_A  = 1/10;             % /           Rapporto area gola / inizio ugello aggiunto

sigma = 0.7;                %              Correction factor for property variations across the boundary layer (all'uscita dal nozzle)

gg = 32.2;                  % ft/s^2       Accelerazione gravitazionale


% OSS: R_d, f_aw, sigma sono presi da pagina 111 del liquid rocket design

% CALCOLI:

T_wg = f_wg * T_c;          % R     Hot-gas-side local chamber-wall temperature
T_aw = f_aw * T_c;          % R     Adiabatic wall temperature of the gas


mu = 46.6e-10 * mm^0.5 * T_c^0.6;      % viscosità
Pr = 4*g / ( 9*g -5 );

hg = (0.026/((R_t*2)^0.2)) * (mu^0.2 * cp / ...
    (Pr^0.6)) * (p_c * gg / c_star)^0.8 * (((R_t*2) / R)^0.1 ) ...
* (rapp_A)^0.9 * sigma;             % gas side heat transfer coefficient [btu/in^2*sec*F]

h_gc = 1 / ( (1/hg) + R_d );        % btu/in^2*sec*F     overall gas-side thermal conductance 

q = h_gc * (T_aw - T_wg);           %[btu/in^2*sec]     Heat flux or heat transferred across the stagnant 
                                    %                   gas film per unit surface area per unit time 

q_real = q * 1055.06/ 0.00064516;   % Heat flux [J/m^2*s]

%% h_g calculation h_g = a*b*c

a = (0.026/((R_t*2)^0.2)) * (mu^0.2 * cp / (Pr^0.6)) * (p_c * gg / c_star)^0.8 * (((R_t*2) / R)^0.1 );
a_vec = a * ones(1,12);

rapp_A_vect = 1./linspace(10,16,12);
b_vec = rapp_A_vect.^0.9;

exp     = linspace(10,16,12);

exp_to_interp  = [10 11 12 13 14 15 16];
r_d_vec_to_int = [2000 2025 2050 2060 2070 2080 2090]; % values from exp 10,11,12,13,14,15,16.

poly = polyfit(exp_to_interp, r_d_vec_to_int, 2);
xx   =  linspace(10,16,100);
yy   = polyval(poly, xx);
% figure;
% plot(xx, yy);
% hold on 
% plot(exp_to_interp, r_d_vec_to_int, 'x');
% grid on;
% grid minor;
r_d_vec   = polyval(poly, exp);

c_vec    =  0.7*ones(1,12); %sigma_vec 

h_g_vec  = a_vec.*b_vec.*c_vec;
 
h_gc_vec = (1./(h_g_vec) + r_d_vec).^-1;

q_conv   = (T_aw - T_wg) * h_gc_vec * 1055.06/ 0.00064516;

%%


T_co  = 1588.25; % °R - temperatura del cooling all'ingresso. 
c_pvc = 0.655;   % btu/(lb F) - specific heat at costant pressure
eta_c = 0.50;    % .25 - .65 valori tipici

rapp = (T_aw - T_wg)/(T_aw - T_co);
esp  = (h_gc_vec./(c_pvc*eta_c));
rapp_2 = 1/rapp;
G_c  = esp./log(rapp_2) *703.0696; % kg/(m^2 sec)
area = 19.9877; % [m^2] internal nozzle area from 10:1 to 16:1
G_total = G_c * area;