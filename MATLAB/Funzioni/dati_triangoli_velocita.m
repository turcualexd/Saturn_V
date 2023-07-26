clear, clc, close all;

c1 = 1.1721e3;                      % m/s
U = 274.1844;                       % m/s
a1 = acos(4*U / c1);                % rad
kn = 0.89;
gamma = 1.128179;
R = 294.4417;                       % J/kg*K
T = 888.38;                         % K

c4 = 0.5 * sqrt(gamma * R * T);     % m/s
a4 = pi/2;

Dh = 7.9290e5;
dh = 0.02 * Dh;

save dati_triangoli_velocita