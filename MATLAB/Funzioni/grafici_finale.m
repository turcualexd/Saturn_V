%% CALCOLO DELLO SCAMBIO TERMICO NELL'UGELLO AGGIUNTIVO
% scambio convettivo, liquid rocket design da pagina 110

clear
close all
clc


% valori nel punto di aggiunta dell'estensione dell'ugello

T_c = 6429*0.975^2; % R            temperatura in camera di combustione con fattore correttivo
p_c = 1124.85;    % psia         

g   = 1.2439;     % /            Gamma dei gas combusti 1:10
mm  = 22.26;       % lb/mol       Massa molare 
cp  = 1.1225;     % btu/lb*F     Specific heat at constant pressure a 1:10
c_star = 5929.2;  % ft/s         

%T_g = 3012.64;    % R            Temperatura gas combusti in 1:10
R_10 = 57.58;     % in           Raggio nel punto di 1:10

f_wg   = 0.8;     %              T_wg/Tc 
%R_d    = 2090;       % in^2*s*F/btu Thermal resistance caused by the solid deposition

R = ((666.75 + 169.80)/2)/25.4;         % in           mean Radius of curvature of nozzle contour (media interno/esterno) (parte nozzle)

    % /           Rapporto area gola / inizio ugello aggiunto

   %              Correction factor for property variations across the boundary layer (all'uscita dal nozzle)

gg = 32.2;        % ft/s^2       Accelerazione gravitazionale


% OSS: R_d, f_aw, sigma sono presi da pagina 111 del liquid rocket design

% CALCOLI:

mu = 46.6e-10 * (mm^0.5) * (T_c^0.6);      % viscosità
Pr_1 = (4*g) / ( 9*g - 5 );

% hg = (0.026/((R_10*2)^0.2)) * (mu^0.2 * cp / ...
%     (Pr^0.6)) * (p_c * gg / c_star)^0.8 * (((R_10*2) / R)^0.1 ) ...
% * (rapp_A)^0.9 * sigma; % gas side heat transfer coefficient [btu/in^2*sec*F]
% 
% h_gc = 1 / ( (1/hg) + R_d ); % btu/in^2*sec*F     overall gas-side thermal conductance 
% 
% q = h_gc * (T_aw - T_wg); %Heat flux or heat transferred across the stagnant 
% % gas film per unit surface area per unit time [btu/in^2*sec]
% 
% q_real = q * 1055.06/ 0.00064516;    % Heat flux [J/m^2*s]

r = Pr_1^0.33;




%% h_g calculation h_g = a*b*c
r_t = 17.5; %inch

a = (0.026/((r_t*2)^0.2)) * ((mu^0.2) * cp / (Pr_1^0.6)) * ((p_c * gg / c_star)^0.8) * ( ((r_t*2) / R)^0.1 );
a_vec = a*ones(100,1);
xx = linspace(0, 112.98, 100)'; %inch - discretizzazione lunghezza ugello sull'asse di simmetria
cc = 0.8*0.445*39.37/tan(deg2rad(15));
rapp = (xx/cc + 1).^-2;
b_vec = rapp.^0.9;

% c_vec    =  [0.95, 0.82, 0.79, 0.78, 0.77,0.76, 0.75, 0.72,0.71 ,0.7]; %sigma_vec 
% poly     =  polyfit(1:10, c_vec, 2);
% y_poly   =  polyval(poly, linspace(1,10,100));
% h_g1_10  =  a_vec.*b_vec.*y_poly;


 
%q_conv   = (T_aw - T_wg) * h_g1_10* 1055.06/ 0.00064516;

%% mach number interpolation (da Ae/At 1 a 10)
A_at  = linspace(1,10,10);
Mach = [1, 2.018, 2.341, 2.552, 2.712, 2.842, 2.951, 3.046, 3.130, 3.205];

figure
grid on;
grid minor;
plot(A_at, Mach, '.k', 'MarkerSize', 8);
xlabel('A_{at}');
ylabel('Mach')
hold on

for i = 1:9
    poly = polyfit(A_at, Mach, i);
    poly_Mach = polyval(poly, A_at);
    [m(i), k(i)] = max(abs(poly_Mach - Mach));
end

poly = polyfit(A_at, Mach, 9);
Aree = rapp.^-1;
poly_Mach = (polyval(poly, Aree));
plot(Aree, poly_Mach, 'green');
legend('data', 'interpolation');

poly_rd = R_d(Aree);
poly_rd = poly_rd';

%%

n = 100; %numero di stazioni tra gola e rapporto 10:1
%T_wg = 0.8*6429*ones(n,1);
T_wg = 5000*ones(n,1);
M = poly_Mach;
T_c = T_c*ones(100,1);
T_ci = zeros(100,1);
T_co = zeros(100,1);
T_wc = zeros(100,1);
T_aw = T_c.*(1 + (r* ((g-1)/2) *(poly_Mach.^2))) ./ (1 + (((g-1)/2) *(poly_Mach.^2)));
T_cool = zeros(n,1);
T_iniziale_coolant = (527 + 600)/2;
T_cool(1) = T_iniziale_coolant; %Rankine -
m_dot_cool = 1145.219881; %lb/s

c_spec_rp1 = 0.5; %btu / lb F 
%R = r_t*sqrt(linspace(1,10,100));
dx = (xx(end)/n);
T_ci(1) = T_iniziale_coolant;
lambda_wall = 3.86*1e-4; %dati design of LRE systems Huzel et al. pag 119 pdfk_cool = 1.78*1e-6; % dati design of LRE ..
mu_cool  = 4.16*1e-5; % dati deisgn  of LRE 
k_cool = 1.78*1e-6;
Pr = (mu_cool*c_spec_rp1/k_cool);

n_tubes_pre  = 178;
n_tubes_post = 356;
diam_tubes_ext   = zeros(100,1);

L_3 = 0.80 * ((sqrt(3) - 1) * 0.4445)/ (tan(deg2rad(15))) * 39.37; %inches (factor 39.37  convert from m to inches)


for k = 1:100

    if xx(k) < L_3

        diam_tubes_ext(k) = 2*pi*(r_t * sqrt(Aree(k))) / n_tubes_pre;
    else
        diam_tubes_ext(k) = 2*pi*(r_t * sqrt(Aree(k))) / n_tubes_post;
    
    end

end

t = 0.0319; %spessore tubi in inches
d_tubes_int = diam_tubes_ext - 2*t;
area_tubes_int = (pi/4) * d_tubes_int.^2;
G = m_dot_cool./area_tubes_int;

% questo modello RPA (ciclo iterativo) non considera il deposito carbonioso del flusso di gas
while k > 0.05
    T_wg(1)
    T_wg_i = T_wg; %per check convergenza
   
    
    sigma_vec = (((0.5*T_wg./T_c) .* (1 + (M.^2) * 0.2439/2) + 1/2).^(-0.68)).*(1 + (M.^2) * 0.2439/2).^(-0.12);
    alfa_t = a_vec.*b_vec.*sigma_vec;
    alfa_t = (alfa_t.^-1 + poly_rd).^-1;
    q_gas_side = alfa_t.*(T_aw - T_wg)
    Radius = r_t*(rapp.^(-1/2));
   
    dT_cool = (2*pi*flip(Radius).*flip(q_gas_side)*dx)./(m_dot_cool*c_spec_rp1);
   
    for i = 1:n-1
        T_ci(i+1) = T_ci(i) + dT_cool(i);
        T_co(i+1) = T_ci(i);
    end
    T_co(1) = T_ci(2);
    T_cool = (T_ci + T_co)/2;
    T_wc = T_wg - q_gas_side*t/lambda_wall;
    h_cool_side = 0.029 * c_spec_rp1 * ((mu_cool^(0.2)) / ((Pr)^(2/3))) .* ( (G.^0.8) ./ (d_tubes_int.^0.2) ) .* (T_cool./T_wc).^0.55;
    %STEP DI ITERAZIONE -> RIDEFINIRE T_WC E T_WG;
    T_wc = flip(T_cool) + q_gas_side./h_cool_side;
    T_wg = (q_gas_side*t/lambda_wall) + T_wc;
    T_wg_i1 = T_wg;
    k = (abs(T_wg_i1 - T_wg_i))./T_wg_i;
    k = min(k);

end

figure;
plot(xx, T_wg,'Color', 'r');
xlabel('asse nozzle [inch]');
ylabel('Temp. wall lato gas T[R]')
hold on;
%plot(xx, T_aw*(5/9),'Color', 'b')
plot(xx, T_wc,'Color', 'g')
xlabel('asse nozzle [inch]');
ylabel('Temp. wall lato coolant T[R]')
plot(xx, flip(T_cool))
xlabel('asse nozzle [inch]');
ylabel('Temp. coolant T[R]')
grid on;
figure;
plot(xx, q_gas_side);
xlabel('asse nozzle [inch]');
ylabel('flusso di calore q[btu/inch^2 s]')

grid on;

%% Pressure drop 
l_nozzle= 3.1898; % nozzle arc length 
dp = zeros(n,1);
rho_rp1 = 0.0293; %lb/inch^3

for i=1:n
   
    if i < 35
        num = 178;
    else 
        num = 356;
    end
    PI = pi*d_tubes_int(i)*num;
    A  = num* pi * d_tubes_int(i)^2 / 4;
    w  = m_dot_cool/(rho_rp1*A);
    Re = rho_rp1 * w * d_tubes_int(i) / mu_cool;
    lambda = 64/Re;
    d_eq = 4*A/PI;
    dp(i) = lambda*(dx/d_eq)*rho_rp1*w^2 / 2;

end

dP = sum(dp)/14.54;


%%


 %T_co  = 1588.25; % °R - temperatura del cooling all'ingresso. 
 %c_pvc = 0.655;   % btu/(lb F) - specific heat at costant pressure
 %eta_c = 0.50;    % .25 - .65 valori tipici
% 
 %rapp = (T_aw - T_wg)/(T_aw - T_co);
% esp  = (h_gc_vec./(c_pvc*eta_c));
 %rapp_2 = 1/rapp;
% G_c  = esp./log(rapp_2) *703.0696; % kg/(m^2 sec)
% area = 19.9877; % [m^2] internal nozzle area from 10:1 to 16:1
 %G_total = G_c * area;
 
 
 %figure;
% plot(xx, q_conv, 'o');
% grid on;
% m_dot = 519.463;
 
% R = [R flip(R)];
% c_spec_rp1 = 2093.399; %j/kgK
% dx = (xx(end)-xx(1))/100;
% q_conv = [q_conv(1:end-1), q_conv];
 
 %for i = 2:199
 %    T(1) = 298.15;
 %    T(i) = T(i-1) + 2*pi*R(i)*dx*(q_conv(i))/(m_dot*c_spec_rp1);
 %end
 %xx = [xx(1:end-1) ,flip(xx)];
 %figure;
 %plot(xx, T, 'o');
 %legend("Cooling Temperature")
 %xlabel("Axis Lenght [m]");
 %ylabel("Cooling RP-1 temperature [K]");
 %grid on;
% grid minor;
 
% q_conv = q_conv(1:100)/(1055.06/ 0.00064516);
% Twg_vec = T_aw*ones(1,100) - q_conv./h_g1_10;
 
% figure;
 %plot(xx(1:100), Twg_vec, 'o');
