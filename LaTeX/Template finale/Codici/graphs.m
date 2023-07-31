clear
close all
clc

temperature = [866.48; 894.26; 922.039; 949.817; 977.59; 1005.37; 1033.15; 1060.93; 1088.71; 1116.48; 1144.26; 1172.04; 1199.82];
cp = [2658.62; 2675.36; 2692.11; 2704.67; 2713.05; 2725.61; 2733.98; 2742.35; 2750.73; 2759.10; 2763.29; 2767.47; 2771.66];
gamma = [1.097; 1.100; 1.106; 1.111; 1.115; 1.119; 1.124; 1.128; 1.132; 1.137; 1.140; 1.144; 1.148];
O_F = [0.308; 0.320; 0.337; 0.354; 0.372; 0.390; 0.408; 0.425; 0.443; 0.460; 0.478; 0.497; 0.516];

%% CP
figure
grid on;
grid minor;
plot(temperature, cp, '.k', 'MarkerSize', 8);
xlabel('Temperature [K]');
ylabel('cp [J/(K kg)]')
hold on

for i = 1:13
    poly = polyfit(temperature, cp, i);
    poly_cp = polyval(poly, temperature);
    [m(i), k(i)] = max(abs(poly_cp - cp));
end
poly = polyfit(temperature, cp, 11);
poly_cp = polyval(poly, temperature);
plot(temperature, poly_cp, 'red');
legend('data', 'interpolation');

%% GAMMA
figure
grid on;
grid minor;
plot(temperature, gamma, '.k', 'MarkerSize', 8);
xlabel('Temperature [K]');
ylabel('\gamma')
hold on

for i = 1:13
    poly = polyfit(temperature, gamma, i);
    poly_gamma = polyval(poly, temperature);
    [m(i), k(i)] = max(abs(poly_gamma - gamma));
end
poly = polyfit(temperature, gamma, 11);
poly_gamma = polyval(poly, temperature);
plot(temperature, poly_gamma, 'blue');
legend('data', 'interpolation');

%% O/F
figure
grid on;
grid minor;
plot(temperature, O_F, '.k', 'MarkerSize', 8);
xlabel('Temperature [K]');
ylabel('O/F')
hold on

for i = 1:13
    poly = polyfit(temperature, O_F, i);
    poly_OF = polyval(poly, temperature);
    [m(i), k(i)] = max(abs(poly_OF - O_F));
end
poly = polyfit(temperature, O_F, 4);
poly_OF = polyval(poly, temperature);
plot(temperature, poly_OF, 'green');
legend('data', 'interpolation');

%% F1
cp_F1 = polyval(polyfit(temperature, cp, 11), 1062);
gamma_F1 = polyval(polyfit(temperature, gamma, 11), 1062);
OF_F1 = polyval(polyfit(temperature, O_F, 4), 1062);