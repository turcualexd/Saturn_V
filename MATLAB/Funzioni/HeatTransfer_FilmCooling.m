%% CALCOLO DELLO SCAMBIO TERMICO NELL'UGELLO AGGIUNTIVO
% scambio convettivo, liquid rocket design da pagina 110

clear
clc
close all 

% valori nel punto di aggiunta dell'estensione dell'ugello

T_c_teorica = 6429.6;            % R            temperatura in camera di combustione
T_c = T_c_teorica*0.975^2;         % R            fattore correttivo di c* pari a 0.975 
p_c = 1124.85;                   % psia         

g   = 1.1754;                    % /            Gamma dei gas combusti 
mm  = 22.26;                     % lb/mol       Massa molare gas combusti (RPA)
cp  = 1.1225;                    % btu/lb*F     Specific heat at constant pressure a 1:10
c_star = 5929.2;                 % ft/s         

R_t = 17.5;                      % in           Raggio in gola

f_wg   = 0.8;                    % T_wg/Tc 
R_d = 2001.8;                    % in^2*s*F/btu Thermal resistance caused by the solid deposition

R = (27.38+6.97)/2;              % in           Radius of curvature of nozzle contour (parte nozzle)

rapp_A  = 1/10;                  % /           Rapporto area gola / inizio ugello aggiunto

sigma = 0.7;                     %              Correction factor for property variations across the boundary layer (all'uscita dal nozzle)

gg = 32.2;                       % ft/s^2       Accelerazione gravitazionale


% OSS: R_d, f_aw, sigma sono presi da pagina 111 del liquid rocket design

% CALCOLI:
Pr = 4*g / ( 9*g -5 );          % considerata costante 

T_wg = f_wg * T_c;              % R     Hot-gas-side local chamber-wall temperature

% calcolo di T_aw = T_c * f_aw         dove f_aw è il turbolent boundary layer recovery factor [val medio 0.9-0.98]

A_at  = linspace(1,10,10);
Mach = [1, 2.018, 2.341, 2.552, 2.712, 2.842, 2.951, 3.046, 3.130, 3.205];

for i = 1:9
    poly = polyfit(A_at, Mach, i);
    poly_Mach = polyval(poly, A_at);
    [m(i), k(i)] = max(abs(poly_Mach - Mach));
end

poly = polyfit(A_at, Mach, 9);
Aree = linspace(10,16,50);
poly_Mach = (polyval(poly, Aree));

r = Pr^0.33;                           % local recovery factor

T_aw = T_c.*(1 + (r* ((g-1)/2) *(poly_Mach.^2))) ./ (1 + (((g-1)/2) *(poly_Mach.^2)));

mu = 46.6e-10 * mm^0.5 * T_c^0.6;       % viscosità

hg = (0.026/((R_t*2)^0.2)) * (mu^0.2 * cp / (Pr^0.6)) * ...
    (p_c * gg / c_star)^0.8 * (((R_t*2) / R)^0.1 ) * ...
    (rapp_A)^0.9 * sigma;               % btu/in^2*sec*F     gas side heat transfer coefficient [btu/in^2*sec*F]

h_gc = 1 / ( (1/hg) + R_d );            % btu/in^2*sec*F     overall gas-side thermal conductance 

q = h_gc * (T_aw - T_wg);               %[btu/in^2*sec] Heat flux or  
                                        % heat transferred across the stagnant gas film per unit surface area per unit time 

q_real = q * 1055.06/ 0.00064516;        % [J/m^2*s]          Heat flux 




