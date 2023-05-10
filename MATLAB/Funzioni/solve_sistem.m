

function F = solve_sistem(x)

C_1   = 1.1721e+03;
U     = 244;
a1    = 0.585685543457151;
kn    = 0.89;
gamma = 1.128179;
R     = 294.4417; %J/kgK
C_4   = 0.4 * (sqrt(gamma*R*888.38));

% x(1) --> alfa_2
% x(2) --> C_2 

F(1)  = (x(2)*cos(x(1) + pi) - U)^2 + (x(2)*sin(x(1) + pi))^2 - (kn^2)*( ((C_1*cos(a1) - U))^2 + (C_1*sin(a1))^2 );

F(2)  = U^2 + C_4^2 - kn^2 * ( (kn*x(2)*cos(x(1)) - U)^2  + (-kn*x(2)*sin(x(1)))^2 );


