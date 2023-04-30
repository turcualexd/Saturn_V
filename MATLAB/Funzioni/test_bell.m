clear
close all
clc

%% Nozzle data:
eps = 16;
Rt = 0.4625;

bell_percentage = 0.8;
theta_n_deg = 33;
theta_e_deg = 7;

[x, y, Ln_cone, Ln, Nx, Ny, Qx, Qy, Ex, Ey] = bell_nozzle(eps, Rt, bell_percentage, theta_n_deg, theta_e_deg);

%% 2D Plot
figure
grid on;
grid minor;
hold on;
plot(x, y);
plot(Nx, Ny, 'o');
plot(Qx, Qy, '.k');
plot(Ex, Ey, 'o');
legend('Bell Nozzle 2D', 'Inflection point', 'Cross point', 'Exit point');

%% Plot 3D
figure
[X, Y, Z] = cylinder(y, 100);
surf(X, Y, Z)

%%
N = length(x) - 1;

theta = linspace(0, 2*pi, 50);

X = zeros(length(theta), N+1);
Y = zeros(length(theta), N+1);
Z = zeros(length(theta), N+1);

for i = 1:length(theta)
    for j = 1:N+1
        X(i,j) = x(j);
        Y(i,j) = y(j)*cos(theta(i));
        Z(i,j) = y(j)*sin(theta(i));
    end
end

figure
surf(X,Y,Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
