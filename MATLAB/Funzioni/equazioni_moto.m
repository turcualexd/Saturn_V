clear
close all
clc

% m0 massa iniziale del sistema
m0 = 0;

% mf massa totale del sistema alla fine della fase propulsa
mf = 0;

% mp massa di propellente necessaria allo svolgimento del sistema
mp = 0;

% MR rapporto di massa del sistema
MR = mf / m0;

% csi_p frazione di massa del propellente
csi_p = 1 - MR;

% m_t variazione di massa nel tempo 
t = 0;       % istante di tempo generico
tb = 0;      %inserire tempo di combustione dello stadio 
m_t = m0 .* (1 - (1-MR) .* (t./tb));

%% Attenzione si usano relazioni differenziali

% c velocità efficace del razzo
c = 0;
% dv = (csi_p ./ (1-csi_p.*(t/tb))).*(c/tb)dt

% Per una traiettoria bi-dimensionale
% A è la superficie dell'intero vettore di lancio che contribuisce a creare drag
% du/dt = ((c*csi_p)/tb)/(1-(csi_p/tb)) - g*sin(theta) - (CD*1/2*pho*u.^2*A/m0)/(1-(csi_p*t/tb)
% sin(theta) = 1 nel caso di traiettoria verticale

% u_p velocità alla fine della combustione del propellente
% u0 è la velocità iniziale (nel nostro caso è nulla)
% B deve essere integrato usando metodi numerici
% B = integrale tra 0 e tb di (1/2*pho*u^2)/(1-csi_p*t/tb)
g = 0;
u_p = -c .* log(1-csi_p) - g*tb - (B*CD*A)/m0 + u0;  % bisogna vedere bene quale valore di g usare

%% Considerazioni su struttura del lanciatore (attenzione perchè le formule non considerano resistenza e gravità in questo caso)

% mp / m0 = 1 - exp(-deltav / (g0*IS))

% deltaV per un singolo stadio
g0 = 9.81;
IS = 0;    % impulso specifico relativo allo stadio considerato
delta_v = g0 .* IS .* log(m0 ./ mf);
