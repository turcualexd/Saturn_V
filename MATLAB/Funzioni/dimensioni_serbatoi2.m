function [lunghezza_rp1, volume_rp1, lunghezza_lox, volume_lox] = dimensioni_serbatoi2(prop_rp1, rho_rp1, prop_lox, rho_lox, diametro, pressione_rp1, pressione_lox)
% Calcola la lunghezza e il volume dei serbatoi del primo stadio del Saturn V
% in base alla quantità e alla densità dei due propellenti, al diametro dello stadio
% e alle pressioni iniziali che si hanno a quota del mare.

% Calcola il volume totale dei propellenti
volume_totale = prop_rp1 / rho_rp1 + prop_lox / rho_lox;

% Calcola il volume di ciascun serbatoio assumendo che abbiano la stessa lunghezza
volume_singolo = volume_totale / 2;

% Calcola il raggio dei serbatoi
raggio = diametro / 2;

% Calcola la lunghezza e il volume del serbatoio RP-1
pressione_critica_rp1 = pressione_rp1 * 4.4; % Conversione da PSI a kPa
lunghezza_rp1 = volume_singolo * pressione_critica_rp1 / (pi * raggio^2 * rho_rp1);
volume_rp1 = lunghezza_rp1 * pi * raggio^2;

% Calcola la lunghezza e il volume del serbatoio LOX
pressione_critica_lox = pressione_lox * 6.9; % Conversione da PSI a kPa
lunghezza_lox = volume_singolo * pressione_critica_lox / (pi * raggio^2 * rho_lox);
volume_lox = lunghezza_lox * pi * raggio^2;
end
