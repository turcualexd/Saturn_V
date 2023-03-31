clear, clc, close all;

%% Parametri e dati iniziali
dt = 0.1;                       % time increment
tb = 161;                       % combustion time
t = 0:dt:tb;                    % time vector
k = length(t);                  % number of iterations

% Costanti e dati a h=0
mu = 3.986005e14;               % gravitational parameter for Earth
R_E =  6373249;                 % Earth radius
T_sl = 33850000;                % thrust sea level
T_vac = 38850000;               % thrust vacuum
I_sp = 304;                     % specific impulse vacuum

g(1) = 9.80665;                 % acceleration of gravity
m(1) = 2898941;                 % total mass
m(k) = 242186;                  % mass dry
m_p(1) = m(1) - m(k);           % mass propellant

[rho(1),a(1),T(1),p(1),nu(1)] = atmos(0);   % standard atmosphere parameters at h=0

m_dot = T_vac / (g(1) * I_sp);  % propellant mass flow rate
A_e = (T_vac - T_sl) / p(1);    % area of nozzle exit


%% Definizione delle funzioni
T_f = @(p) T_vac - A_e .* p;    % thrust
T(1) = T_f(p(1));               % T at h=0

