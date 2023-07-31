clear; close all; clc;

%F-1 engine specs
d_e = 3.53;
A_e = pi * d_e^2 /4 ;
m_dot_o = 1788.97 + 22.23;
m_dot_f = 742.09  + 53.52;
m_dot = m_dot_f + m_dot_o;
R     = 8314;            %J/kg K
g     = 1.174;
M_mm  = 22.2152;
p_e   = 0.0423; %MPa
T_c   = 3500; %combustion temp

v_e   = @(p_c) sqrt(2*R*(g/(g-1))*(T_c/M_mm)*(1-(p_e./p_c).^((g-1)/g)) ); %m/s
gg = 9.81;

pp   = linspace(2,45,1000);
vv   = v_e(pp);

%system spec
%fuel and lox pumps
dP_ox =      11.045;
k_ox = dP_ox/(7.7566);
rho_o = 1145;
dP_f  =      12.893;
k_f = dP_f/(7.7566);
rho_f = 810;
eta_o_pump = 0.746;
eta_f_pump = 0.726;

%turbine
eta_t = 0.605;
T_in  = 1061;
%from CEA hP problem (gas gen)
M_mm_gg = 19.247;
 P_te = 0.3998;
gamma = 1.128179;
c_p_gas =  2742.2380;

P_c_gg = 0.85 * pp; %best practice 
eta_tt = 1 - (P_te./P_c_gg).^((gamma-1)/(gamma));

dP_oo = pp*k_ox;
dP_ff = pp*k_f;

%power balance
Pw_lox = (m_dot_o*dP_oo*1e6)/(eta_o_pump*rho_o);
Pw_rp1 = (m_dot_f*dP_ff*1e6)/(eta_f_pump*rho_f);
req_pw = Pw_lox + Pw_rp1;

m_dot_gg = req_pw./(eta_t*c_p_gas*T_in*eta_tt);

v_e_gg   = @(p_c) sqrt( 2*R*(gamma/(gamma-1))*(T_in/M_mm_gg)*(1-(P_te./P_c_gg).^((gamma-1)/gamma)));

I_gg = v_e_gg(P_c_gg)/gg;
m_chamb = m_dot*ones(1, 1000) - m_dot_gg; 
I_spec_see = @(v, m_ch) (v + (A_e./m_ch)*1e6*(p_e - 0.101325))/gg;
I_ss       = I_spec_see(vv, m_chamb);
I_tot = (m_dot_gg/m_dot).*I_gg + (1 - m_dot_gg/m_dot).*I_ss;


% specific impulse see level
figure;
title('Total system specific impulse at s.l.')
plot(pp, I_tot,'LineWidth',1.5);
xlabel('chamber pressure [MPa]');
ylabel('Specific Impulse at see [s]');
hold on
plot(7.765*ones(1000,1), linspace(250, 320, 1000), '--r');
plot(pp, I_ss, 'LineWidth',1.5)
legend('I_{s,oa}', 'Duty Point', 'I_{s,tc}');