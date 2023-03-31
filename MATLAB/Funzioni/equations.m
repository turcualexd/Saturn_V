clear
close all
clc

g_0= 9.80665; 
dt = 0.1;            % time increment
tb = 161;            % combustion time

m_t = 0;
m_t0 = 2898941;      % total mass
m_d = 242186;        % mass dry
m_p = m_t0 - m_d;    % mass propellant

T_sl = 33850000;     % thrust sea level
T_vac = 38850000;    % thrust vacuum
I_sp = 304;          % specific impulse vacuum

m_dot = T_vac / (g_0 * I_sp);  % propellant mass flow rate

p_a0 = 101325;
A_e = (T_vac - T_sl) ./ p_a0;  % area of nozzle exit

T_0 = T_vac - A_e .* p_a0;     % thrust

% definire angoli di pitch theta da tabella 

mu = 3.986005e14;   % gravitational parameter for Earth
R_E =  6373249;

% definire h quota del vettore di lancio durante tutte le varie fasi

for t = 0:dt:tb

    m_p = m_p - m_dot .* dt;
    m_t = m_d + m_p;

    % T = T_vac - A_e .* p; !! Attenzione p dipende dalla quota h(i)
    % R = R_E + h(i);
    % g = mu ./ R;

    % definire accelerazioni e velocità
    % 16 equazioni da risolvere 


end