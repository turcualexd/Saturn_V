clear
close all
clc

% Dati stadio:
r2_r1 = 1;
cm2_cm1 = 1;
D0s_D0h = 1.25;
n = 25000;
alpha1 = deg2rad(75);
u = 350;
eta_r = 0.92;
eta_s = 0.88;
f = 0.03;

% Fluido di lavoro:
T00 = 0;
p00 = 0;
k = 1.4;
R_f = 287;
cp = 1004;

psi = [1e-1:0.1:1e3]';
R = [0];

eta_tts = [];

for i = 1:1
    for j = 1:size(psi,1)
        kis = 1./sqrt(psi(j));
        C1 = sqrt(eta_s)/kis .* (sqrt(1-R(i)));
        C1t = C1.*sin(alpha1);
        C2m = C1.*sin(alpha1);
        C1m = C2m;
        W2 = sqrt(((eta_r*(1+f).*R(i))/kis.^2) + (eta_s/(kis.^2)).*(1-R) - 2*sqrt(eta_s)/(kis).*sqrt(1-R(i))*sin(alpha1) + (r2_r1).^2);
        W2t = sqrt(W2.^2 - C2m.^2);
        C2t = (W2t - r2_r1);

        C2 = sqrt(C2m.^2 + C2t.^2);
        lambda = 2*(C1t - r2_r1*C2t);

        phie = 1;

        eta = lambda/(psi(j) - phie*C2.^2);
        eta_tts = [eta_tts; eta];

    end
end

%%
figure
hold on;
grid on;
grid minor;
plot(psi, eta_tts);
