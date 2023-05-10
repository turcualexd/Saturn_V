clear, clc, close all;

c1 = 1.1721e3;                      % m/s
U = 244.1844;                       % m/s
a1 = acos(4*U / c1);                % rad
kn = 0.89;
gamma = 1.128179;
R = 294.4417;                       % J/kg*K
T = 888.38;                         % K

c4 = 0.65 * sqrt(gamma * R * T);     % m/s
a4 = pi/2;

x0 = [600, pi/6];

x = fsolve(@sist, x0);

c2 = x(1);      % m/s
a2 = x(2);      % rad

c3 = kn * c2;   % m/s
a3 = a2;        % rad