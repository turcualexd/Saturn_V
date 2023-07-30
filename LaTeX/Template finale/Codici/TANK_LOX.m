clear;clc;
%% LOX

rhoLOX = 1141; % kg/m^3
g = 9.81;

Vtp = 44716; %ft^3 volume propellente
Vt_p = 77269485;%in^3
a = 198; %inch %raggio della sezione cilindrica del serbatoio
b = 120; %inch altezza cupola elissoidale
lc = 528; %inch %altezza parte cilindrica
NPSH = 720; %inch %vedi pompe

%% volumi

V_e = (2*pi*a^2*b)/3;  %inch^3
V_c = pi*a^2*lc;    %inch^3
V_he = 4*53568; % volume occupato da elio %inch^3
V_tot = 2*V_e + V_c - V_he;  %inch^3

V_u = V_tot - Vt_p;  %inch^3
Vu_1 = V_e - V_u; %prova che l'ullage è minore del volume dell'ellissoide, quindi bu meno di 
syms b_2
H_u = solve( V_u == (2*pi*a^2*b_2)/3, b_2); % altezza ullage inch  

%% pressioni

H_in = lc + 2*b - H_u; %inch
H_prop = H_in /39.37; %m

p_u = 23 ; %pressione ullage +-2.5 psia

pi_n = rhoLOX*g*H_prop; %Pascal
%p_i = pi_n /6895 % psia

p_i = 27.92; %psia
p_t = p_i + p_u; %psia
%pt = 54.00755736; %psi


%% materiale

rho = 0.103; %lb/in^3
s_r = 69000; %psi
s_y = 57000; %psi
E = 10.6e6; %psi
v = 0.33; %poisson
% sicurezza
S_r = s_r/1.3;
S_y = s_y/1.25;



%% spessori

k = a/b; %elipse ratio
K = 0.8; % stress factor
R = a*k;

t_e = (p_t*a*(K+0.5*k))/(2*S_y); %spessore ellissoide %inch
t_c = (p_t*a)/(S_y); %spessore cilindro %inch

%%

t_e_m = t_e*25.4; %mm 
t_c_m = t_c*25.4; %mm

%% aree

e = (sqrt(a^2-b^2))/a; % eccentricità

A_e = a^2 + (pi*b^2*log((1+e)/(1-e)))/(2*e); %area ellissoidale

A_c = 2*pi*a*lc;

A_tot = 2*A_e + A_c; %inch^2

%% peso


E1 = 2*k+(1/sqrt(k^2-1))*log((k+sqrt(k^2-1))/(k-sqrt(k^2-1))); %design factor
W_e = (pi*a^2*t_e*E1*rho)/(2*k); 

W_c = 2*pi*a*lc*t_c*rho; 

PesoLOX = 2*W_e + W_c; %lb