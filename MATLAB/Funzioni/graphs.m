clear
close all
clc

temperature = [866.48; 894.26; 922.039; 949.817; 977.59; 1005.37; 1033.15; 1060.93; 1088.71; 1116.48; 1144.26; 1172.04; 1199.82];
cp = [2658.62; 2675.36; 2692.11; 2704.67; 2713.05; 2725.61; 2733.98; 2742.35; 2750.73; 2759.10; 2763.29; 2767.47; 2771.66];
gamma = [1.097; 1.100; 1.106; 1.111; 1.115; 1.119; 1.124; 1.128; 1.132; 1.137; 1.140; 1.144; 1.148];
O_F = [0.308; 0.320; 0.337; 0.354; 0.372; 0.390; 0.408; 0.425; 0.443; 0.460; 0.478; 0.497; 0.516];

figure;
grid on;
grid minor;
tx = 800:0.001:1250;
cpq1 = interp1(temperature, cp, tx,'spline');
plot(temperature, cp,'o', tx, cpq1, ':.');

figure;
grid on;
grid minor;
tx = 800:0.001:1250;
gammaq1 = interp1(temperature, gamma, tx,'spline');
plot(temperature, gamma,'o', tx, gammaq1, ':.');

figure;
grid on;
grid minor;
tx = 800:0.001:1250;
O_Fq1 = interp1(temperature, O_F, tx);
plot(temperature, O_F,'o', tx, O_Fq1, ':.');