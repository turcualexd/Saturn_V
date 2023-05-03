clear
close all
clc

%% Nozzle 1:10 data:
eps = 10;
Rt = 0.4625;

bell_percentage = 0.8;
[theta_n_deg, theta_e_deg] = bellvalues(bell_percentage, eps);

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

%% 3D Plot
n = length(x);
theta_3d = linspace(0, 2*pi, n);
X = repmat(x, [1,n]);
Y = y.*cos(theta_3d);
Z = y.*sin(theta_3d);

figure
surf(X, Y, Z);
axis equal
box on
xlabel('X');
ylabel('Y');
zlabel('Z');

%% Nozzle 1:16 data:
eps = 16;
Rt = 0.4625;

bell_percentage = 0.8;
[theta_n_deg, theta_e_deg] = bellvalues(bell_percentage, eps);

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

%% 3D Plot
n = length(x);
theta_3d = linspace(0, 2*pi, n);
X = repmat(x, [1,n]);
Y = y.*cos(theta_3d);
Z = y.*sin(theta_3d);

figure
surf(X, Y, Z);
axis equal
box on
xlabel('X');
ylabel('Y');
zlabel('Z');