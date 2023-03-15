clear; clc; close all;

p_e = 0.6619*1e5; %Pa
p_c = 7.756600*1e6;
T_c = 3572; %K
Mm  = 22.372;
gamma_GC = 1.1484;
m_prop = 2601.81; %kg/s

u_e = exit_v(gamma_GC, Mm, T_c, p_e, p_c);
p_array = standard_atmosphere (0:200:6000);
pe = p_e * ones(31,1);
A_e = pi*(3.7/2)^2;
T = Thrust(m_prop, u_e, pe', p_array,A_e);
z = linspace(0,6000, 31);
plot(z, T)