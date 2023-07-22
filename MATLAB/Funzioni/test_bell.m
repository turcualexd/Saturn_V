clear
close all
clc

%% Nozzle 1:10 data:
eps_10 = 10;
% throat diameter 35 inches = 0.889 m
Rt_10 = 0.4445;

bell_percentage_eps_10 = 0.8;
[theta_n_deg_eps_10, theta_e_deg_eps_10] = bellvalues(bell_percentage_eps_10, eps_10);

[x_eps_10, y_eps_10, Ln_cone_eps_10, Ln_eps_10, Nx_eps_10, Ny_eps_10, Qx_eps_10, Qy_eps_10, Ex_eps_10, Ey_eps_10] = bell_nozzle(eps_10, Rt_10, bell_percentage_eps_10, theta_n_deg_eps_10, theta_e_deg_eps_10);

si = 0;
for i=2:125
    si = si + sqrt((x_eps_10(i)-x_eps_10(i-1)).^2 + (y_eps_10(i)-y_eps_10(i-1)).^2);
end

%% 2D Plot
figure
grid on;
grid minor;
hold on;
plot(x_eps_10, y_eps_10);
plot(Nx_eps_10, Ny_eps_10, 'o');
plot(Qx_eps_10, Qy_eps_10, '.k');
plot(Ex_eps_10, Ey_eps_10, 'o');
legend('Bell Nozzle 2D', 'Inflection point', 'Cross point', 'Exit point');

%% 3D Plot
n = length(x_eps_10);
theta_3d = linspace(0, 2*pi, n);
X_eps_10 = repmat(x_eps_10, [1,n]);
Y_eps_10 = y_eps_10.*cos(theta_3d);
Z_eps_10 = y_eps_10.*sin(theta_3d);

figure
surf(X_eps_10, Y_eps_10, Z_eps_10);
axis equal
box on
xlabel('X');
ylabel('Y');
zlabel('Z');

%% Nozzle 1:16 data:
eps = 16;
% throat diameter 35 inches = 0.889 m
Rt = 0.4445;

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
