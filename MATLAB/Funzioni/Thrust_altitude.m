clear; clc; close all;

p_e = 0.6619*1e5; %Pa
p_c = 7.756600*1e6;
T_c = 3572; %K
Mm  = 22.372;
gamma_GC = 1.1484;
m_prop = 2601.81; %kg/s

u_e = exit_v(gamma_GC, Mm, T_c, p_e, p_c);
