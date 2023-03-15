function [P_z] = standard_atmosphere (altitude)

if (altitude >= 0 || altitude < 11000)
    % valori di riferimento
    g = 9.81;      %[m/s]
    P_0 = 101325;  %[Pa]
    T_0 = 288.15;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    % gradiente di temperatura
    a = 0.0065;     %[K/m]

    P_z = P_0 .* (1 - ((a*altitude)/T_0)).^(g/(R*a));
end

if (altitude >= 11000 || altitude <= 20000)
    % valori di riferimento
    z_tropo = 11000;
    g = 9.81;      %[m/s]
    P_tropo = 22626;  %[Pa]
    T_tropo = 216.65;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    P_z = P_tropo .* exp(-(g/(R.*T_tropo)).*(altitude - z_tropo));
end

if (altitude >= 20000 || altitude < 32000)
    % valori di riferimento
    g = 9.81;      %[m/s]
    P_0 = 101325;  %[Pa]
    T_0 = 288.15;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    % gradiente di temperatura
    a = 0.001;     %[K/m]

    P_z = P_0 .* (1 - ((a*altitude)/T_0)).^(g/(R*a));
end

if (altitude >= 32000 || altitude < 47000)
    % valori di riferimento
    g = 9.81;      %[m/s]
    P_0 = 101325;  %[Pa]
    T_0 = 288.15;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    % gradiente di temperatura
    a = 0.0028;     %[K/m]

    P_z = P_0 .* (1 - ((a*altitude)/T_0)).^(g/(R*a));
end

if (altitude >= 47000 || altitude <= 51000)
    % valori di riferimento
    z_tropo = 47000;
    g = 9.81;      %[m/s]
    a = 0.0028;
    P_tropo = P_0 .* (1 - ((a*47000)/T_0)).^(g/(R*a));;  %[Pa]
    T_tropo = 270.65;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    P_z = P_tropo .* exp(-(g/(R.*T_tropo)).*(altitude - z_tropo));
end

if (altitude >= 51000 || altitude <= 61000)
    % valori di riferimento
    g = 9.81;      %[m/s]
    P_0 = 101325;  %[Pa]
    T_0 = 288.15;  %[K]

    % costante specifica dell'aria
    R = 287.1;      %[J/(kg*K)]

    % gradiente di temperatura
    a = -0.0028;     %[K/m]

    P_z = P_0 .* (1 - ((a*altitude)/T_0)).^(g/(R*a));
end


