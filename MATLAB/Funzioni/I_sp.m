clear; close all; clc;

%F-1 engine specs
A_e = 10.75210;
m_dot_o = 1844;
m_dot_f = 777;
m_dot = 1844 + 777;
R     = 8314;            %J/kg K
g     = 1.174;
M_mm  = 22.2152;
p_e   = 0.049; %MPa
T_c   = 3500; %combustion temp
v_e   = @(p_c) sqrt(2*R*(g/(g-1))*(T_c/M_mm)*(1-(p_e./p_c).^((g-1)/g)) ); %m/s
gg = 9.81;
I_spec_see = @(v) (v + (A_e/m_dot)*1e6*(p_e - 0.101325))/gg;
I_spec_vac = @(v) (v + (A_e/m_dot)*p_e*1e6)/gg;

pp   = linspace(2,10,1000); %from 0.2 to 10Mpa. Pressure at chamber 7.57 MPa
vv   = v_e(pp); 
I_ss = I_spec_see(vv);
I_ss_vac = I_spec_vac(vv);



%%

%hold on
%set = 0.77*1e7*ones(100,1);
%yy  = linspace(100, 400,100);
%plot(set, yy, 'r--');
%GG CYCLE

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
% gamma_s = 1.1378;
% der = -1.12404;
% gamma= - gamma_s*der;
 P_te = 0.3998;
% c_p_gas = (gamma/(gamma - 1)) * R/M_mm_gg;
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

I_tot = (m_dot_gg/m_dot).*I_gg + (1 - m_dot_gg/m_dot).*I_ss;
I_tot_vac = (m_dot_gg/m_dot).*I_gg + (1 - m_dot_gg/m_dot).*I_ss_vac;

%required mass flow

figure;
title('mass flux in GG')
plot(pp, m_dot_gg);
xlabel("chamber pressure [MPa]");
ylabel("mass flux gg [kg/s]");
hold on
plot(7.765*ones(100,1), linspace(0, 90, 100))
legend('', 'Duty Point');

% specific impulse see level

figure;
title('Total system specific impulse at s.l.')
plot(pp, I_tot);
xlabel('chamber pressure [MPa]');
ylabel('Specific Impulse at see [s]');
hold on
plot(7.765*ones(1000,1), linspace(220, 350, 1000))
plot(pp, I_ss)
legend('', 'Duty Point', 'theoretical $I_{spec}$');

%gg chamber pressure

figure;
title('chamber pressure in GG')
plot(pp, P_c_gg);
xlabel("main chamber pressure [MPa]");
ylabel("gg chamber pressure [MPa]");
hold on
plot(7.765*ones(100,1), linspace(0, 10, 100))
legend('', 'Duty Point');

% vacum

figure;
title('Total system specific impulse vacuum')
plot(pp, I_tot_vac);
xlabel('chamber pressure [MPa]');
ylabel('Specific Impulse vacuum [s]');
hold on
plot(7.765*ones(1000,1), linspace(200, 350, 1000))
legend('', 'Duty Point');



