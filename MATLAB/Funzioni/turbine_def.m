clear; close all; clc;

% dati iniziali (interpolazione dati AIAA + manuali F-1)
c_p       = 2742.2380 ;      % c_p gas generator from interpolation
T_gg      = 1062;  % 
T_turb_e  = 888.38;
p_in      = 63.43; % bar
eps       = 16.4;
p_e       = p_in/eps;
o_f       = 0.416; %
g         = 1.128179;
R         = 294.4417;
mdot      = 77.92;
omega     = 574.7;
U_over_c0 = 0.225;

% ipotesi nostre sy rendimenti (AIAA)
eta_ovrll = 0.6;    % rendimento totale turbina (termodinamic + meccanico)
noz_ar    = 9.7;    % aspect ratio nozzle
k_n       = 0.96;   % perdita spouting velocity
eta_n     = k_n^2;  % rendimento ugello
eps_nt    = 0.97;   %area coefficient throat
eps_ne    = 0.95;   %area cofficient exit
k_blade   = 0.89;   % rotor and stator blade velocity coefficient
eps_blade = 0.95;   % stator and rotor blade exit area coefficient


% struttura turbina 

% nozzle inlet - nozzle out --> rotor --> stator --> rotor

%% Inlet nozzle -> suppongo che il 94% del salto di dh avvenga nei nozzle

dh_tot     = c_p*T_gg * (1 - (1/eps)^((g-1)/g)); % salto totale turbina.
dh_nz      = (1 - 0.06)*dh_tot;
dh_rot1    = dh_tot*0.02;
dh_stat    = dh_rot1;
p_e_nzz    = p_in*(1 - dh_nz/(c_p*T_gg))^(g/(g-1));
C_zero     = sqrt(2*dh_nz);
C_1        = k_n * C_zero; 
U          = C_zero*U_over_c0;

U_vec       = U*[1, 0];
q_n_re      = ((1 - k_n^2)*C_1^2)/(k_n^2);
T_e_nzz_id  = T_gg - dh_nz/c_p;
T_e_nzz_re  = T_e_nzz_id + q_n_re/c_p;

%relazione per efficienza blade massima
alfa_1     = acos(4*U/C_1);
alfa_1_deg = alfa_1*(180/pi);

beta_1     = atan(C_1*sin(alfa_1)/(C_1*cos(alfa_1) - U));
beta_1_deg = rad2deg(beta_1);


C_1_vec    = C_1 * [cos(-alfa_1), sin(-alfa_1)];
V_1_vec    = C_1_vec - U_vec;
V_1        = norm(V_1_vec);

%% solving system 
% ipotesi --> mach uscita 0.4, alfa2 = alfa3, c3 = kc2, v2= kv1 e v4 = kv3
% risolvendo il sistema troviamo alfa 2 e C2

x0     = [600, pi/4];
x_sol  = fsolve(@sist, x0); % troviamo rispettivamente C_2 e alfa2
C_2    = x_sol(1);
alfa_2 = x_sol(2);

%% first rotor exit - stator  inlet


T_e_rot1_id = T_e_nzz_re - dh_rot1/c_p;
q_rot1_re   = (1 - k_blade^2)*(V_1^2 / 2) + (1 - k_n^2)*dh_rot1;
T_e_rot1_re = T_e_rot1_id + q_rot1_re/c_p;
C_2_vec     = [C_2*cos(pi + alfa_2), C_2*sin(pi + alfa_2)];
V_2_vec     = C_2_vec - U_vec;
p_e_rot1    = p_e_nzz*(1 - dh_nz/(c_p*T_e_nzz_re))^(g/(g-1));
beta_2      = acos(V_2_vec(1)/norm(V_2_vec));

%% stator outlet - 2nd rotor inlet

alfa_3      = alfa_2;
C_3         = k_blade*C_2;
C_3_vec     = C_3*[cos(-alfa_3), sin(-alfa_3)];
V_3_vec     = C_3_vec - U_vec;
q_bs        = (1 - k_blade^2)*(C_2^2 / 2) + (1 - k_n^2)*dh_stat;
T_e_stat_id = T_e_rot1_re - dh_stat/c_p;
T_e_stat_re = T_e_stat_id + q_bs/c_p;
p_e_stat    = p_e_rot1*(1 - dh_stat/(c_p*T_e_rot1_re))^(g/(g-1));


%% 2nd rotor outlet

alfa_4      = pi/2;
C_4         = 0.4*sqrt(g*R*T_turb_e);
C_4_vec     = C_4 * [cos(-alfa_4), sin(-alfa_4)];
V_4_vec     = C_4_vec - U_vec;
q_rot2_re   = (1 - k_blade^2)*(norm(V_3_vec)^2 / 2) + (1 - eta_n)*dh_rot1;
T_e_rot2_id = T_e_stat_re - dh_rot1/c_p;
T_e_rot2_re = T_e_rot2_id + q_rot2_re/c_p;

%% diagrammi

figure;

pos1 = [0,0,C_1_vec(1), C_1_vec(2)];
quiver2(pos1,'b');

hold on;

pos2 = [0,0, V_1_vec(1), V_1_vec(2)];
quiver2(pos2,'r');

pos3 = [V_1_vec(1), V_1_vec(2), U, 0];
quiver2(pos3,'#38DB4A');

grid on
%grid minor
plot(linspace(-200, 1200, 10000), C_1_vec(2)*ones(10000,1), 'LineStyle', '--', color="k");

pos4 = [V_1_vec(1), V_1_vec(2), C_2_vec(1), C_2_vec(2)];
quiver2(pos4,'b');

pos5 = [V_1_vec(1), V_1_vec(2), V_2_vec(1), V_2_vec(2)];
quiver2(pos5,'r');

pos6 = [V_1_vec(1) + V_2_vec(1), V_1_vec(2) + V_2_vec(2), U, 0];
quiver2(pos6,'#38DB4A');

pos7 = [V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), C_3_vec(1), C_3_vec(2)];
quiver2(pos7,'b');

pos8 = [V_1_vec(1) + C_2_vec(1), V_1_vec(2) + C_2_vec(2), V_3_vec(1), V_3_vec(2)];
quiver2(pos8,'r');

pos9 = [V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), U, 0];
quiver2(pos9,'#38DB4A');

pos10 = [V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), C_4_vec(1), C_4_vec(2)];
quiver2(pos10,'b');

pos11 = [V_1_vec(1) + C_2_vec(1) + V_3_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2), V_4_vec(1), V_4_vec(2)];
quiver2(pos11,'r');

pos12 = [V_1_vec(1) + C_2_vec(1) + V_3_vec(1) + V_4_vec(1), V_1_vec(2) + C_2_vec(2) + V_3_vec(2) + V_4_vec(2), U, 0];
quiver2(pos12,'#38DB4A');

plot(linspace(-200, 1200, 10000),(C_1_vec(2) + C_2_vec(2))*ones(10000,1),'LineStyle', '--', color="k");
plot(linspace(-200, 1200, 10000), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2))*ones(10000,1), 'LineStyle', '--', color="k");
plot(linspace(-200, 1200, 10000), (C_1_vec(2) + C_2_vec(2) + C_3_vec(2) + C_4_vec(2))*ones(10000,1), 'LineStyle', '--', color="k");
plot(0,0,'k+')

text(650,-150,'C_1','FontSize',18,'FontWeight','bold')
text(650,-650,'C_2','FontSize',18,'FontWeight','bold')
text(600,-950,'C_3','FontSize',18,'FontWeight','bold')
text(550,-1300,'C_4','FontSize',18,'FontWeight','bold')
text(300,-250,'V_1','FontSize',18,'FontWeight','bold')
text(300,-550,'V_2','FontSize',18,'FontWeight','bold')
text(310,-1000,'V_3','FontSize',18,'FontWeight','bold')
text(230,-1250,'V_4','FontSize',18,'FontWeight','bold')
text(900,-500,'U_1','FontSize',18,'FontWeight','bold')
text(130,-900,'U_2','FontSize',18,'FontWeight','bold')
text(510,-1090,'U_3','FontSize',18,'FontWeight','bold')
text(380,-1320,'U_4','FontSize',18,'FontWeight','bold')

%% nozzle dimensioning

A_nozz_throat = mdot/(eps_nt*p_in*1e5*sqrt( (g*(2/(g+1))^((g+1)/(g-1)))/(R*T_gg) ) ); %m^2 di area di gola totale
H_nozz_throat = 1.55*0.0254; % inch --> to m (dal manuale si ha altezza di 1.55 inch
b_nozz        = (H_nozz_throat/noz_ar); %larghezza nozzle --> si ipotizza un Aspect Ratio di 9.7
z_n           = A_nozz_throat/(H_nozz_throat*b_nozz); %stima numero nozzle (totali veri 61, qui calcolati 69 --> diminuire aspect ratio)

%% efficiency 

eta_nb_prac = (U*(norm(C_1_vec)*cos(alfa_1) + norm(C_2_vec)*cos(alfa_2) + ...
    norm(C_3_vec)*cos(alfa_3) + norm(C_4_vec)*cos(alfa_4))) / dh_tot;

eta_m       = eta_ovrll/eta_nb_prac;


function quiver2(pos, col)
    a = annotation('arrow', 'HeadStyle', 'plain', 'HeadLength', 5, 'HeadWidth', 5, 'LineWidth', 1, 'Color', col);
    set(a, 'parent', gca);
    set(a, 'position', pos);
end