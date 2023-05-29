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
R_t = 17.5;       % in           Raggio di gola

f_aw   = 0.923;   % /            Turbolent boundary layer recovery factor in gola [val medio 0.9-0.98]
R_d = 1125;       % in^2*s*F/btu Thermal resistance caused by the solid deposition in gola

R = 6.96;         % in           Radius of curvature of nozzle contour (parte nozzle)
rapp_A  = 1;      % /            Rapporto aree in gola

gg = 32.2;        % ft/s^2       Accelerazione gravitazionale


% calcoli preliminari all'altezza della gola:

T_wg = 1188;            % R     Hot-gas-side local chamber-wall temperature (gola) -> dato sperimentale (p111 modern engine)
T_aw = f_aw * T_c;      % R     Adiabatic wall temperature of the gas

sigma = 1.42;           % /     Correction factor for property variations across the boundary layer (in gola): dal grafico p110 liquid rocket design, considerando T_wg/T_c = 0.1848 e interpolando (approx) 


mu_gg= 46.6e-10 * mm^0.5 * T_g^0.6;     % viscosità gas combusti
Pr = 4*g / ( 9*g -5 );                  % considero costante 

hg = (0.026/((R_t*2)^0.2)) * (mu_gg^0.2 * cp / (Pr^0.6)) * (p_c * gg / c_star)^0.8 * (((R_t*2) / R)^0.1 ) * (rapp_A)^0.9 * sigma; % gas side heat transfer coefficient [btu/in^2*sec*F]
h_gc = 1 / ( (1/hg) + R_d );            % btu/in^2*sec*F     overall gas-side thermal conductance 

q = h_gc * (T_aw - T_wg);               % btu/in^2*sec       Heat flux or heat transferred across the stagnant gas film per unit surface area per unit time 


%% REGENERATIVE

%DATI

% dati lega X750
t = 0.0298;          % in  *************    Thickness (per far venire corretto il numero di tubi
k_lega = 3.19e-4;    % btu/in^2 *s * F/in   Thermal conductivity of chamber wall 

E = 28e6;            % psi                  Modulus of elasticity of tube wall material
a = 8e-6;            % in/in * F            Thermal expansion coefficient of tube wall material
v = 0.35;            % 7                    Coefficiente di Poisson

% dati RP1 tubi in gola
T_co   = 600;        % R    ********          Coolant bulk temperature in gola nei tubi up ( temperatura accettabile per i baffles e piu alta della temperatura di stoccaggio)
                     %                        ben sotto la temperatura critica di ebollizione del kerosene

k_fuel = 1.78e-6;    % btu/in^2*s*(F/in)      Coolant thermal conductivity 
C1 = 0.0214;         % costante propria del RP1 per il calcolo del numero di nusselt

% RP1 in condizioni T_co
mu   =  4.16e-5;     % lb/in*s                Coolant viscosity at bulk temperature 600 R
mu_w =  0.416e-5;    % lb/in*s    ********    Coolant viscosity at coolant sidewall temperature (valore preso dall'esempio, con Twc = 1000, questo è 994)         

w_f = 70/100 * 1754.02;                 % lb/s         mass flow rate RP1 che passa nei tubi
ro =  0.0293;                           % lb/in^3      densità RP1 (valore di stivaggio, approssimazione nella trattazione pompa: no delta ro no delta T

% calcoli: dimensionamento tubi altezza gola considerando tubo up (condizione peggiore) 

T_wc = T_wg - (q*t / k_lega);           % R             coolant side wall temperature
h_c  = q / (T_wc - T_co);               % Btu/in^2sF    coefficiente di scambio termico per consentire flusso termico per il salto di temperatura dalla t parete lato fluido e la t media fluido   

% Nusselt = C1*Re^0.*  Pr* 0.4 * rapp_mu ^ 0.14
% sostituendo : hc*d/k_f = C1 * ( ro*V_co*d / mu) * a * b  con d diametro interno tubi incognita
% cerco coefficienti aa bb, V_co 

aa = (mu * cp / k_fuel) ^ 0.4; 
bb = (mu / mu_w) ^ 0.14;

% sostituisco l'espressione di V_co funzione di N e d nell'equazione dei
% Nusselt per trovare N = N(d)

syms d N_tubi            % d diametro interno, N = numero tubi

V_co = (w_f/ro) / (N_tubi/2 * (pi*d^2)/4);  % in/s         coolant velocity

Coeff_N = ( 0.0214 * ( 8*w_f/(mu*pi) ) ^ 0.8 * (aa * bb * k_fuel / h_c) ) ^ (1/0.8); 
% dove N è il numero di tubi calcolato come Coeff_N * diametro?(-9/4)

% Sostituisco N = pi*(R_t*2 + 0.8*(diametro_interno+0.04))/(diametro_interno+0.04) 
% nell'espressione per il calcolo del numero di tubi e trovo d

L = 0.8 * 0.04 + 2*R_t;     

% risolvendo (Coeff_N * d^(-9/4) == pi* (0.8*d + L) / (d +  0.04)) si
% ottiene la soluzione in modo grafico (wolframalpha) oppure tramite solve

f(d) = Coeff_N * d^(-9/4) - pi* (0.8*d + L) / (d +  0.04);
[sol, par, cond] = solve( f(d) == 0, d, 'Real', true, 'ReturnConditions',true);

assume(cond)
restr = sol > 0;
sol_par = solve(restr,par);

d_interno = double(subs(sol,par,sol_par));                % in

N_tubi = round(Coeff_N * d_interno^(-9/4));         
diametro_interno = d_interno * 2.54;                      % cm

V_co = (w_f/ro) / (N_tubi/2 * (pi*d_interno^2)/4);        % in/s


%% Dimensionamento a stress 

% DATI
p_co = 1500;         % psia     p_RP1 intermedia tra p_uscita dalla pompa fuel (1860 psia) 
                     %          e quella dell'injector-manifold (1222 psia)

T_amb = -9.67;       % F        temperatura ambiente alla massima quota raggiunta dal primo stadio
T_c_F = 5970;        % F        temperatura in camera di combustione in F
% CALCOLI

p_g = p_c * (2/(g+1))^(g/(g-1));     % psia   pressione gas combusti in gola

% Tensione di trazione tangenziale (t) e longitudinale (l)

S_t = ((p_co - p_g)*(d_interno/2)/t) * (E*a*q*t) / (2*(1-v)*k_lega);    % lb/in^2

%***********
S_l = E * a * ( T_c - T_amb);   % lb/in^2 

% S_l è stato determinato al massimo salto di temperatura (condizione piu
% critica) e in queste condizioni questo stress non può essere maggiore
% dello 0.9 dello stress critico. da questa osservazione si puo calcolare
% lo stress critico:

S_c = S_l/0.9; % stress critico per instabilità anelastica longitudinale

%***********

% caduta di pressione per unità di lunghezza all'altezza della gola

Re = ro * V_co * d_interno / mu;   
lambda = 0.0032 + 0.221/(Re^0.237);        % fattore di attrito per turbolent flows
Delta_p = lambda * ro * V_co^2 / (d_interno * (gg*12)); %lb/in^3  

h_c
