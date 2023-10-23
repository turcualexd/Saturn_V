clear, clc, close all;

h = linspace(0, 62000, 100000);
[~,~,~,p] = atmos(h);

figure
plot(h, p)
hold on

% for i = 1:40
%     poly = polyfit(h, log(p), i);
%     poly_h = polyval(poly, h);
%     [m(i), k(i)] = max(abs(exp(poly_h) - p));
% end
% il polinomio di grado 13 è la scelta più ponderata

poly = polyfit(h, log(p), 13);
poly_h = polyval(poly, h);
plot(h, exp(poly_h))
