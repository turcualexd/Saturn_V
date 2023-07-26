clear; clc; close all;
t_burn = 168;
p_e = 0.6619*1e5; %Pa
p_c = 7.756600*1e6;
T_c = 3572; %K
Mm  = 22.372;
gamma_GC = 1.1484;
m_prop = 2601.81; %kg/s

u_e = exit_v(gamma_GC, Mm, T_c, p_e, p_c);
z = linspace(0,60000, 30000);
[~,~,~,P] = atmos(z);

pe = p_e * ones(30000,1);
A_e = pi*(3.7/2)^2;
T = Thrust(m_prop, u_e, pe', P, A_e);
plot(z, P);
figure
plot(z, T, "Color", "r");
xlabel("quota [m]"); ylabel("Spinta [N]");
grid on;

I_sp = T./(9.81*m_prop);
% figure
% plot(z, I_sp)

M_i = 2.278059*1e6;
t_burn;
tt = linspace(0, t_burn, 1000);
M = @(t) M_i - 5*m_prop * t;
MM = M(tt);
figure;
plot(tt, MM)
xlabel("time [s]");
ylabel("massa [kg]");
grid on;
MM(end)

