function [lunghezza, volume] = dimensioni_serbatoi(prop1, rho1, prop2, rho2, diametro)
% Calcola la lunghezza e il volume dei serbatoi del primo stadio del Saturn V
% in base alla quantità e alla densità dei due propellenti e al diametro dello stadio.

% Calcola il volume totale dei propellenti
volume_totale = prop1 / rho1 + prop2 / rho2;

% Calcola il volume di ciascun serbatoio assumendo che abbiano la stessa lunghezza
volume_singolo = volume_totale / 2;

% Calcola il raggio dei serbatoi
raggio = diametro / 2;

% Calcola la lunghezza di ciascun serbatoio
lunghezza = volume_singolo / (pi * raggio^2);

% Calcola il volume di ciascun serbatoio
volume = lunghezza * pi * raggio^2;
end
