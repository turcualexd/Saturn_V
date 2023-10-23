% dimensioning casing of turbopumps

clear; close all; clc

Q      = 25200;
k_v    = 0.35;
H      = 3.099831932773109e+03;
c3     = k_v*sqrt(2*32.2*H);
alfa_2 = 0.152980330947031;

a = @(theta) (theta*Q)./(3.12*360*c3);

tt = linspace(0,180,90);
a_tt = a(tt);

a_v    = 2*a_tt(end);
alfa_v = alfa_2;
r_tong = 1.05 * 0.49530/2;
d_in   = sqrt((4*a_v)/pi);
a_in   = a_v;
d_ext  = d_in + 2*10*tan(deg2rad(5));
a_ext  = d_ext^2 * pi/4;

v_in   = Q/(3.12 * a_in);
v_ext  = Q/(3.12 * a_ext);


