function F = sist(x)

load("dati_triangoli_velocita.mat")

f_c_mag = @(c,a) (c * sin(a))^2 + (c * cos(a) - U)^2;
f_c_min = @(c,a) (c * sin(a))^2 + (c * cos(a) + U)^2;

F(1) = f_c_min(x(1), x(2)) - kn^2 * f_c_mag(c1, a1) - 2 * dh;
F(2) = f_c_min(c4, a4) - kn^2 * f_c_mag(kn*x(1), x(2)) - 2 * dh;