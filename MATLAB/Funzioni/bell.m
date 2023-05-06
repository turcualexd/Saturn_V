clear
close all
clc

eps = 16; % expansion ratio of the nozzle
Rt = 0.4625;  % radius of the throat
theta_n = deg2rad(33);
theta_e = deg2rad(7);

Ln_cone = ((sqrt(eps) - 1) * Rt)/ (tan(deg2rad(15)));

Re = sqrt(eps)*Rt;

% 80% bell
Ln = 0.8*((sqrt(eps) - 1) * Rt)/ (tan(deg2rad(15)));

x = [];
y = [];

theta_i = [deg2rad(-135):0.01:deg2rad(-90)];
for i = 1:length(theta_i)
    x_i = 1.5*Rt*cos(theta_i(i));
    y_i = 1.5*Rt*sin(theta_i(i)) + 1.5*Rt + Rt;
    x = [x; x_i];
    y = [y; y_i];
end

theta_exit = [deg2rad(-90):0.01:deg2rad(33-90)];
for i = 1:length(theta_exit)
    x_i = 0.382*Rt*cos(theta_exit(i));
    y_i = 0.382*Rt*sin(theta_exit(i)) + 0.382*Rt + Rt;
    x = [x; x_i];
    y = [y; y_i];
end

Nx = x(end,1);
Ny = y(end,1);

Ex = Ln;
Ey = Re;

m1 = tan(theta_n);
m2 = tan(theta_e);
C1 = Ny - m1*Nx;
C2 = Ey - m2*Ex;

Qx = (C2-C1) / (m1-m2);
Qy = (m1*C2 - m2*C1) / (m1 - m2);

t = [0:0.01:1];
for i = 1:length(t)
    x_b = (1-t(i)).^2 * Nx + 2*(1-t(i))*t(i) * Qx + t(i).^2 * Ex;
    y_b = (1-t(i)).^2 * Ny + 2*(1-t(i))*t(i) * Qy + t(i).^2 * Ey;
    x = [x; x_b];
    y = [y; y_b];
end


figure
grid on;
grid minor;
hold on;
plot(x, y);
plot(Nx, Ny, 'o');
plot(Ex, Ey, 'o');
plot(Qx, Qy, 'o');

%% Plot 3D
figure
[X, Y, Z] = cylinder(y, 100);
surf(X, Y, Z, 'EdgeColor', 'none')

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
surf(X,Y,Z,'EdgeColor','none');
xlabel('x');
ylabel('y');
zlabel('z');