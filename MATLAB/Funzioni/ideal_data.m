function ideal = ideal_data(gamma, MM, T_amb, p_amb,T_c, p_e, A_e, p_c, A_t, u_e)

R = 8314;  % [kJ / kg * K]
g0 = 9.81  % [m/s^2]

% Vandencherkove function
GAMMA = sqrt(gamma.*(2/(gamma+1)).^((gamma+1)/(gamma-1)));

% C_star
c_star = (p_c .* A_t) ./ (GAMMA .* (p_c ./ (sqrt((R / MM).*T_c))) .* A_t);

% CT_ideal
ct_ideal = sqrt((2.*gamma.^2)/(gamma-1).*(2/(gamma+1)).^((gamma+1)/(gamma-1) .* ...
            (1-(p_e/p_c).^((gamma-1)/gamma)))) + ((p_e - p_amb)/p_c).*(A_e/A_t);

% IS_ideal
IS_ideal = (ct_ideal .* c_star) ./ (g0);

% Thrust
T = mp * u_e + (p_e - p_amb) * A_e;

% mass flowrate ideal
m_p_dot = T ./ (IS_ideal .* g0);

% T_efflux
T_e = T_c .* (p_e./ p_c) .^ ((gamma-1)./gamma);

% sound velocity
c_e = sqrt( gamma * (R./MM) *T_e);

% Mach_efflux
M_e = u_e ./ c_e;

% areas
eps = 1 ./ M_e .* [2 ./ (gamma +1) .* (1 + (gamma -1)./2 .* Me_e.^2)].^...
      ((gamma + 1) ./ (2 .* (gamma -1)));
