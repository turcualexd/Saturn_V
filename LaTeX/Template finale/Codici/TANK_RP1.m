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

Ve = (2*pi*a^2*b)/3; 
Vc = pi*a^2*lc;
l_t = lc + 2*b;
d_t = 17;
V_t = 5*pi*(d_t/2)^2*l_t; %volume approssimato occupato dai tubi di lox

Vtot = 2*Ve + Vc - V_t;

Vu = Vtot - Vt_p;
Vu_1 = Ve - Vu; %prova che l'ullage è minore del volume dell'ellissoide, quindi bu meno di b
syms b_2
bu = solve( Vu == (2*pi*a^2*b_2)/3, b_2); % altezza ullage inch  90.9363

%% pressioni
H_in = lc + 2*b - bu; %inch
H = H_in/39.37;
 
pu = 25.2 ; %pressione ullage psia

pi_n = rhoRP1*g*H;
%p_i = pi_n /6895 % psia

p_i = 13.38; %psia
pt = p_i + pu; %psia

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

tk = (K*pt*a)/(S_y); %spessore nodo
tcr = (pt*R)/(2*S_y); %spessore corona

te = (tk+tcr)/2; %spessore equivalente di un ellissoidale (pt*a*(K+0.5*k))/(2*S_y)
tc = (pt*a)/(S_y); %spessore cilindro

%% aree

e = (sqrt(a^2-b^2))/a; % eccentricità

Ae = a^2 + (pi*b^2*log((1+e)/(1-e)))/(2*e); %area ellissoidale

Ac = 2*pi*a*lc;

Atot = 2*Ae + Ac;

%% peso


E1 = 2*k+(1/sqrt(k^2-1))*log((k+sqrt(k^2-1))/(k-sqrt(k^2-1))); %design factor
We = (pi*a^2*te*E1*rho)/(2*k); 

Wc = 2*pi*a*lc*tc*rho; 


%% pressioni critiche dovute ai carichi esterni

Cb = 0.10; %buckling coefficient da 0.05 a 0.10
Pcre = (Cb*2*E*te^2)/(R^2); 

if lc < 4.9*a*sqrt(a/tc)
    Pcrc = 0.807*((E*tc^2)/(lc*a))*((1/(1-v^2))^3*(tc^2/a^2))^(1/4);
else
    Pcrc = (E*tc^3)/(4*(1-v^2)*a^3);
end


%%


tk_m = tk/39.37
tcr_m = tcr/39.37
te_m = te/39.37 
tc_m = tc/39.37





