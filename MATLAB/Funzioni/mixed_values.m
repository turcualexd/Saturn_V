clear
close all
clc

% 1 dati relativi all'ugello 10:1

m1_dot = 2601.81;   %[kg/s]
p1 = 1.004e5;       %[bar]
T1 = 1919.41;       %[K]
MM1 = 23.135;       %[kg/kmol]
gamma1 = 1.2174;
cp1 = 2.0246;
cv1 = cp1/gamma1;

% 2 dati relativi ad exhaust da turbina

m2_dot = 75.75;     %[kg/s]
p2 = 4e5;           %[bar]
T2 = 1056.41;       %[K]
MM2 = 17.58;        %[kg/kmol]
gamma2 = 1.1333;
cp2 = 12.2221;
cv2 = cp2/gamma2;

% mixed
m_mixed = m1_dot + m2_dot;
p_mixed = (p1*m1_dot + p2*m2_dot) / (m1_dot + m2_dot);
T_mixed = (T1*m1_dot + T2*m2_dot) / (m1_dot + m2_dot);
MM_mixed = (MM1*m1_dot + MM2*m2_dot) / (m1_dot + m2_dot);
gamma_mixed = (gamma1*m1_dot + gamma2*m2_dot) / (m1_dot + m2_dot);

cp_mixed = (cp1*m1_dot + cp2*m2_dot) / (m1_dot + m2_dot);
cv_mixed = (cv1*m1_dot + cv2*m2_dot) / (m1_dot + m2_dot);
gamma_mixed_giusto = cp_mixed / cv_mixed;
T_mixed_giusto = (m1_dot*cp1*T1 + m2_dot*cp2*T2) / (m1_dot*cp1+ m2_dot*cp2);
