clear, clc, close all;

%% Initial data and parameters
dt = 0.01;                      % time increment
tb = 161;                       % combustion time
t = 0:dt:tb;                    % time vector
k = length(t);                  % number of iterations

mu = 3.986005e14;               % gravitational parameter for Earth
R_E =  6373249;                 % Earth radius
T_sl = 33850967;                % thrust sea level
diam_e = 3.53;                  % diameter exit noozle
A_e = diam_e^2 * pi / 4 * 5;    % area of nozzle exit
I_sp = 304;                     % specific impulse vacuum
S = 113;                        % drag surface
g0 = 9.80665;                   % acceleration of gravity
m_d = 242186;                   % mass dry
m_i = 2898941;                  % total mass at t=0
mp_i = m_i - m_d;               % propellant mass at t=0

% propellant mass
mp = nan(1,k);
mp(1) = mp_i;

% total mass
m = nan(1,k);
m(1) = m_i;

% height
h = nan(1,k);
h(1) = 0;

% horizontal shift
s = nan(1,k);
s(1) = 0;

% velocity
v = nan(1,k);
v(1) = 0;
v_v = nan(1,k);
v_v(1) = 0;
v_h = nan(1,k);
v_h(1) = 0;

% standard atmosphere
rho = nan(1,k);
c = nan(1,k);
Temp = nan(1,k);
p = nan(1,k);
nu = nan(1,k);
[rho(1),c(1),Temp(1),p(1),nu(1)] = atmos(0);

T_vac = T_sl + p(1) * A_e;      % thrust vacuum
m_dot_5 = T_vac / (g0 * I_sp);  % propellant mass rate 5 motors

% Mach number
Ma = nan(1,k);
Ma(1) = 0;

% acceleration of gravity
g = nan(1,k);
g(1) = mu / (R_E + h(1))^2;

% thrust
T = nan(1,k);
T(1) = T_vac - A_e * p(1);

% drag
D = nan(1,k);
D(1) = 1/2 * rho(1) * v(1)^2 * S * drag_coeff(Ma(1));

% acceleration
a = nan(1,k);
a(1) = -g(1) + (T(1) - D(1)) / m(1);
a_v = nan(1,k);
a_v(1) = a(1);
a_h = nan(1,k);
a_h(1) = 0;

theta = pitch(t);

% flight angle
phi = nan(1,k);
phi(1) = 0;


%% Solution
for i = 2:k
    h(i) = h(i-1) + v_v(i-1)*dt;
    s(i) = s(i-1) + v_h(i-1)*dt;

    v_v(i) = v_v(i-1) + a_v(i-1)*dt;
    v_h(i) = v_h(i-1) + a_h(i-1)*dt;
    v(i) = v(i-1) + a(i-1)*dt;

    phi(i) = atan(v_h(i) / v_v(i));

    [rho(i),c(i),Temp(i),p(i),nu(i)] = atmos(h(i));

    Ma(i) = v(i) / c(i);
    g(i) = mu / (R_E + h(i))^2;
    D(i) = 1/2 * rho(i) * v(i)^2 * S * drag_coeff(Ma(i));

    if t(i) <= 135
        mp(i) = mp(i-1) - m_dot_5 * dt;
        m(i) = m(i-1) - m_dot_5 * dt;
        T(i) = T_vac - A_e .* p(i);
    else
        mp(i) = mp(i-1) - 4/5 * m_dot_5 * dt;
        m(i) = m(i-1) - 4/5 * m_dot_5 * dt;
        T(i) = 4/5 * T_vac - A_e .* p(i);
    end

    a_v(i) = -g(i) + (T(i) * cos(theta(i)) - D(i) * cos(phi(i))) / m(i);
    a_h(i) = (T(i) * sin(theta(i)) - D(i) * sin(phi(i))) / m(i);
    a(i) = sqrt(a_v(i)^2 + a_h(i)^2);

end

%% Plots

% height
figure
plot(t, h)
xlabel("t [s]")
ylabel("h [m]")
title("Quota di volo")
xline(135, '--r')

% velocity
figure
plot(t, v)
xlabel("t [s]")
ylabel("v [m/s]")
title("Velocità di volo")
xline(135, '--r')

% temperature
% figure
% plot(t, Temp)
% xlabel("t [s]")
% ylabel("T [K]")
% title("Temperatura ambientale")

% pressure
% figure
% plot(t, p)
% xlabel("t [s]")
% ylabel("p [Pa]")
% title("Pressione ambientale")

% Mach number
% figure
% plot(t, Ma)
% xlabel("t [s]")
% ylabel("Ma [-]")
% title("Numero di Mach")

% acceleration of gravity
% figure
% plot(t, g)
% xlabel("t [s]")
% ylabel("g [m/s^2]")
% title("Accelerazione di gravità")

% thrust
figure
plot(t, T)
xlabel("t [s]")
ylabel("T [N]")
title("Spinta dei motori")
xline(135, '--r')

% drag
figure
plot(t, D)
xlabel("t [s]")
ylabel("D [N]")
title("Drag aerodinamico")
xline(135, '--r')

% acceleration
figure
plot(t, a)
xlabel("t [s]")
ylabel("a [m/s^2]")
title("Accelerazione del razzo")
xline(135, '--r')

% total mass
figure
plot(t, m)
xlabel("t [s]")
ylabel("m [kg]")
title("Massa totale del razzo")
xline(135, '--r')

% trajectory
figure
plot(s, h, 'r')
xlabel("x [m]")
ylabel("y [m]")
title("Traiettoria di volo")