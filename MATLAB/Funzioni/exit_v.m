function [u_e] = exit_v(gamma, Mmol, TC, Pe, Pc )

Ru = 8314; %KJ/Kg mol

u_e = sqrt((2*gamma)/(gamma - 1)*(Ru/Mmol)*TC*(1 - (Pe/Pc)^(1 - 1/gamma)));

end

