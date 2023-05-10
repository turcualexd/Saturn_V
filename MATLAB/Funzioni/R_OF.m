clear, clc, close all;

conv = 5.38032045589557;
R = [50.4 51.8 53.6 55.4 58.0 59.0 60.7 62.4 64] * conv;
OF = [0.372 0.390 0.408 0.425 0.443 0.460 0.478 0.497 0.516];
n = length(OF);

x = linspace(OF(1), OF(n), 100000);

% for i = 1:40
%     poly = polyfit(h, log(p), i);
%     poly_h = polyval(poly, h);
%     [m(i), k(i)] = max(abs(exp(poly_h) - p));
% end
% il polinomio di grado 13 è la scelta più ponderata

plot(OF, R, 'x')
hold on

poly = polyfit(OF, R, 2);
poly_x = polyval(poly, x);
plot(x, poly_x)

OF_star = 0.416;
R_star = polyval(poly, OF_star)