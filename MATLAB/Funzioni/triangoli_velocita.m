clear, clc, close all;

load("dati_triangoli_velocita.mat")

x0 = [600, pi/4];
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000);
x = fsolve(@sist, x0, options);

c2 = x(1);      % m/s
a2 = x(2);      % rad

c3 = kn * c2;   % m/s
a3 = a2;        % rad

rad2deg(a2)