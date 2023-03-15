function [P_z] = standard_atmosphere (altitude)

% valori di riferimento
g = 9.81;      %[m/s]
P_0 = 101325;  %[Pa]
T_0 = 288.15;  %[K]

% costante specifica dell'aria
R = 287.1;      %[J/(kg*K)]

% gradiente di temperatura
a = 0.0065;     %[K/m]

P_z = P_0 .* (1 - ((a*altitude)/T_0)).^(g/(R*a));

end