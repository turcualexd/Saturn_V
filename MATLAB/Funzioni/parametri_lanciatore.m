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

% M_0 = M_p + M_s + M_u

M_0 = 2898941;
M_u = 659293;
M_s = 824511;
M_p = 2074430;
M_m = 8391.4588;
M_tank = 124799.5412;

% eps_u: indice di massa del carico utile
eps_u = M_u / (M_0 - M_u);

% csi_u: rapporto di massa del carico utile
csi_u = M_u / M_0;

% eps_s: indice strutturale delle masse inerti
eps_s = M_s / (M_0 - M_u);

% MR: mass ratio
MR = eps_s * (1 - csi_u) + csi_u;

% csi_p: frazione di propellente

csi_p = M_p / M_0;

% csi_m: rapporto di massa
csi_m = M_m / M_0;

% rapporto M_tank / M_p
rapporto = M_tank / M_p;
