function F = sist(x)

c1 = 1.1721e3;                      % m/s
U = 244.1844;                       % m/s
a1 = acos(4*U / c1);                % rad
kn = 0.89;
gamma = 1.128179;
R = 294.4417;                       % J/kg*K
T = 888.38;                         % K
c4 = 0.65 * sqrt(gamma * R * T);     % m/s
a4 = pi/2;

Dh = 7.9290e5;
dh = 0.02 * Dh;

f_c_mag = @(c,a) (c * sin(a))^2 + (c * cos(a) - U)^2;
f_c_min = @(c,a) (c * sin(a))^2 + (c * cos(a) + U)^2;

F(1) = f_c_min(x(1), x(2)) - kn^2 * f_c_mag(c1, a1) - 2 * dh;
F(2) = f_c_min(c4, a4) - kn^2 * f_c_mag(kn*x(1), x(2)) - 2 * dh;