clear
close all
clc

%% Formule valide per una singola pompa

% g accelereazione di gravità
g = 9.81;
% Rendimento pompa
eta_p = 0;
% Pressione totale all'uscita dalla pompa in psia
P02 = 0;
% Pressione totale all'ingresso della pompa in psia
P01 = 0;
% Densità fluido ingresso della pompa in lbm/ft^3
pho1 = 0;
% Q portata volumetrica
Q = 0;
% R2 raggio impeller
R2 = 0;
% b2 spessore impeller
b2 = 0;
% omega velocità angolare impeller
omega = 0;
% D diametro impeller
D = 0;

% P1 pressione inlet
P1 = 0;
% Psat pressione saturazione - passaggio di stato liquido vapore
Psat = 0;
% pho densità fluido ingresso inlet
pho = 0;
% Cm coefficiente di momento
Cm = 0;
% R1 raggio inlet
R1 = 0;

%% Design pompa

% H Head rise: parametro usato per definire le performance della pompa
H_feet = (144.*(P02 - PO1) ./ pho1); % in feet da convertire
H_m = H_feet / 3.281;

% delta_h altezza di sollevamento
delta_h = (g .* H_m) ./ eta_p;

% vr2 velocità radiale di uscita del fluido da impeller
vr2 = Q / (2*pi*R2*b2);

% Potenza impeller e quindi gruppo pompa si ricava angolo di deflessione
beta2 = atan((omega.*R2)./(vr2) .* (1 - (H_m*g)./(eta_p.*(omega.*R2).^2)));

% psi Head coefficient
psi = eta_p .* (1 - (vr2./(omega.*R2)).*tan(beta2));

% phi2 Flow coefficient
phi2 = Q ./ (2*pi.*(R2./b2).*(omega.*R2));

%% Parametri adimensionali

% ds specific diameter
ds = (D*(g*H_m).^(1/4))/(Q.^(1/2));

% ns specific speed
ns = (omega*(Q.^(1/2)))/((g*H_m).^(3/4));

prod = ns*ds;

%% Cavitazione

% NPSP net positive suction pressure
NPSP = P1 - Psat;

% NPSH net positive suction head
NPSH = NPSP / (pho*g);

% Thoma parameter sigma
% parametro empirico per cavitazione
sigma = NPSP / (1/2 * pho * Cm.^2);

phi_t = Cm / (omega*R1);

% Ss suction specific speed
Ss = 2.981 / (sigma.^(3/4) .* phi_t);

% Metodo alternativo per cavitazione
tao = NPSP/(1/2*pho*(omega*R1).^2);
theta = atan(phi_t);
Zt = sin(theta) / (1+cos(theta)) * phi_t;

% tao = 3*Zt come verifica;