clear, clc, close all;

%data
mpunto_f = 742.09;      %kg/s
mpunto_ox = 1788.97;    %kg/s
deltap_f = 641;         %KPa 
deltap_ox = 2100;       %KPa
nhole_f = 1404;
nhole_ox = 1428;
ro_f = 810;             %kg/m³
ro_ox = 1141;           %kg/m³
A_f = 0.05484;          %m^2
A_ox = 0.03968;         %m^2


A_1H_f = A_f/nhole_f;
A_1H_ox = A_ox/nhole_ox;

D_1H_f = sqrt((4*A_1H_f)/pi);
D_1H_ox = sqrt((4*A_1H_ox)/pi);

cd_f = mpunto_f/(A_f * sqrt(2*deltap_f*ro_f));
cd_ox = mpunto_ox/(A_ox * sqrt(2*deltap_ox*ro_ox));

%da bernoulli
v_f = cd_f*sqrt(2*deltap_f/ro_f);       %17.1 m/s
v_ox = cd_ox*sqrt(2*deltap_ox/ro_ox);   %40.5 m/s