function [H_m, delta_h, psi, P, ds, ns, prod] = turbopump(eta_p, P01, P02, pho1, Q, m_dot, Dt2, omega)

g = 9.81;

% Design della pompa
% eta_p:   rendimento pompa
% P01:     pressione totale all'ingresso della pompa in psia
% P02:     pressione totale all'uscita dalla pompa in psia
% pho1:    densità fluido ingresso della pompa in lbm/ft^3
% Q:       portata volumetrica
% m_dot:   portata massica
% Dt2:     diametro impeller
% R2:      raggio impeller
R2 = Dt2 / 2;
% b2:      spessore impeller
% beta2:   angolo di deflessione
% omega:   velocità angolare impeller

% Cavitazione 
% P1:      pressione inlet
% Psat:    pressione saturazione - pressione di vapore
% pho:     densità fluido ingresso inlet
% R1:      raggio inlet
% Dt1 = rapportoD * Dt2;
% R1 = Dt1/2;
% theta:   angolo di leading-edge dell'induttore

% Design
% H Head rise: parametro usato per definire le performance della pompa
H_feet = (144.*(P02 - P01) ./ pho1); % in feet da convertire
H_m = H_feet / 3.281;

% delta_h altezza di sollevamento
delta_h = (g .* H_m) ./ eta_p;

% Power
P = m_dot * (g*H_m)/eta_p;

% % Potenza impeller e quindi gruppo pompa si ricava angolo di deflessione
% beta2 = atan((omega.*R2)./(vr2) .* (1 - (P)./(m_dot.*(omega.*R2).^2)));

% psi Head coefficient
% psi = eta_p .* (1 - (vr2./(omega.*R2)).*tan(beta2));
psi = H_m / ((omega*R2).^2 / g);

% Parametri adimensionali
% ds specific diameter
ds = (Dt2*(g*H_m).^(1/4))/(Q.^(1/2));

% ns specific speed
ns = (omega*(Q.^(1/2)))/((g*H_m).^(3/4));

prod = ns*ds;

% prod_controllo = 2 / sqrt(psi);