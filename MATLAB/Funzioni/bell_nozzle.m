function [x, y, Ln_cone, Ln, Nx, Ny, Qx, Qy, Ex, Ey] = bell_nozzle(eps, Rt, bell_percentage, theta_n_deg, theta_e_deg)

theta_n = deg2rad(theta_n_deg);
theta_e = deg2rad(theta_e_deg);

% Length of a standard 15-degree half angle conical nozzle
Ln_cone = ((sqrt(eps) - 1) * Rt)/ (tan(deg2rad(15)));

% Radius of the nozzle exit
Re = sqrt(eps)*Rt;

% Nozzle length for a __% bell
Ln = bell_percentage * ((sqrt(eps) - 1) * Rt)/ (tan(deg2rad(15)));

x = [];
y = [];

% Equation for the throat:

% Entrant section
theta_i = [deg2rad(-135):0.01:deg2rad(-90)];
for i = 1:length(theta_i)
    x_i = 1.5*Rt*cos(theta_i(i));
    y_i = 1.5*Rt*sin(theta_i(i)) + 1.5*Rt + Rt;
    x = [x; x_i];
    y = [y; y_i];
end

% Exit section
theta_exit = [deg2rad(-90):0.01:deg2rad(theta_n-90)];
for i = 1:length(theta_exit)
    x_i = 0.382*Rt*cos(theta_exit(i));
    y_i = 0.382*Rt*sin(theta_exit(i)) + 0.382*Rt + Rt;
    x = [x; x_i];
    y = [y; y_i];
end

% Points for the mathematical costruction of the bell:

% N: inflection point
Nx = x(end,1);
Ny = y(end,1);

% E: exit point of the nozzle
Ex = Ln;
Ey = Re;

% Q: cross point
m1 = tan(theta_n);
m2 = tan(theta_e);
C1 = Ny - m1*Nx;
C2 = Ey - m2*Ex;
Qx = (C2-C1) / (m1-m2);
Qy = (m1*C2 - m2*C1) / (m1 - m2);

% Equation for the bell: quadratic Bezier curve:

t = [0:0.01:1];
for i = 1:length(t)
    x_b = (1-t(i)).^2 * Nx + 2*(1-t(i))*t(i) * Qx + t(i).^2 * Ex;
    y_b = (1-t(i)).^2 * Ny + 2*(1-t(i))*t(i) * Qy + t(i).^2 * Ey;
    x = [x; x_b];
    y = [y; y_b];
end
