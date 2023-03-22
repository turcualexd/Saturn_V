clc; close all; clear;

FILE = fopen("SI_C.txt", "r");
formatSpec = '%f';

A = fscanf(FILE, formatSpec,[3,inf])';

fclose(FILE);

AeAt = A(:, 1);
Pe   = A(:, 2);
Te   = A(:, 3);

figure;
plot(AeAt, Pe, 'o');
grid on;
grid minor;
xlabel("Rapporto di Ae/At");
ylabel("Pressione [bar]");
