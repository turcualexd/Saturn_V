clear; clc;
%% RP-1

rhoRP1 = 810; %kg/m^3; 
g = 9.81;

Vtp = 27939; %ft^3 volume propellente
Vt_p = 48278996.2;%in^3
a = 198; %inch
b = 120; %inch
lc = 276; %inch
NPSH = 660; %inch

%% volumi

V_e = (2*pi*a^2*b)/3;  %inch^3
V_c = pi*a^2*lc;  %inch^3
l_tubi = lc + 2*b;  %inch
d_tubi = 17;  %inch
V_tubi = 5*pi*(d_tubi/2)^2*l_tubi; %volume approssimato occupato dai tubi di lox %inch^3

V_tot = 2*V_e + V_c - V_tubi;  %inch^3

V_u = V_tot - Vt_p;  %inch^3
Vu_1 = V_e - V_u; %prova che l'ullage è minore del volume dell'ellissoide, quindi bu meno di b
syms b_2
H_u = solve( V_u == (2*pi*a^2*b_2)/3, b_2); % altezza ullage inch  90.9363

%% pressioni
H_in = lc + 2*b - H_u; %inch
H = H_in/39.37; %m
 
p_u = 25.2 ; %pressione ullage psia

pi_n = rhoRP1*g*H;
%p_i = pi_n /6895 % psia

p_i = 13.38; %psia
p_t = p_i + p_u; %psia

%pt = 40.13; %psi

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

t_e = (p_t*a*(K+0.5*k))/(2*S_y); %spessore elissoide %inch
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

PesoRP1 = 2*W_e + W_c; %lb








