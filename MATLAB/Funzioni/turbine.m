clear; close all; clc;

%interpolazione dati
%..

c_p     = 2742.2380 ;      % c_p gas generator from interpolation
T_gg    = 1062;  % K
p_in    = 63.43; % bar
eps     = 16.4;
p_e     = p_in/eps;
o_f     = 0.416; %
g       = 1.128179;
C_zero  = sqrt(2*c_p*T_gg * (1 - (1/eps)^((g-1)/g) ));
mdot    = 77.92;
omega   = 574.7;
eta     = 0.605;
pwr     = 39.45*1e6; %watt

%pitch line velocity (omega * R_turbina)
U      = 256; %m/s
r_turb = U/omega; %raggio turbina-rotore (torna con manuale nasa)
vel_ratio = U/C_zero;
psi = 0.5*(C_zero/U)^2;


eta_t = @(k) pwr*(md(k)*c_p*T_gg*(1 - (1./k).^((g-1)/g))).^-1;
% kk    = linspace(2,30, 100);
% figure
% plot(kk, eta_t(kk));
% hold on
% plot(16.4*ones(100,1), linspace(0,1,100));

% isoentripic dh

dh_iso_tot = c_p * T_gg * (1 - (1/eps)^((g-1)/g));
