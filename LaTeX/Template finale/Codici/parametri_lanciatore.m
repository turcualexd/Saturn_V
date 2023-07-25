clear
close all
clc

% Dati:
% M_0: massa iniziale del sistema lanciatore             [kg]
% M_u: massa carico utile del lanciatore                 [kg]
%      coincide con la massa degli stadi superiori escluso il primo S-IC
% M_s: massa della struttura del lanciatore              [kg]
% M_p: massa del propellente per il primo stadio S-IC    [kg]
% M_m: massa motore                                      [kg]
% M_tank: massa tank                                     [kg]

% Dati iniziali ricavati da simulazione lancio e manuale F1 Rocketdyne
M_0 = 2898941;         
M_u = 824511;
M_s = 242186;
M_p = 2074430;         % 2656755 - 582325 kg
M_m = 8445.4363;       % 18619 lb - manuale F1 schema pagina 10-
M_tank = 124799.5412;

% eps_u: indice di massa del carico utile
eps_u = M_u / (M_0 - M_u);

% csi_u: rapporto di massa del carico utile
zeta_u = M_u / M_0;

% eps_s: indice strutturale delle masse inerti
eps_s = M_s / (M_0 - M_u);

% MR: mass ratio
MR = (M_0 - M_p) / M_0 ;

% csi_p: frazione di propellente
zeta_p = M_p / M_0;

% csi_m: rapporto di massa motore
zeta_motore = 5*M_m / M_0;

% rapporto M_tank / M_p
rapporto_tank = M_tank / M_p;
