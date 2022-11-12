Boltz = 5.67*10^-8; % W/m^2K^4
Emissivity = .84;
l_sec = 19; %in
l_sec_m = l_sec * 0.0254; %m
d_sec = 6.17; %in
d_sec_m = d_sec * 0.0254; %m
T_amb_F = 110; %F
T_amb = (T_amb_F - 32) * 5/9 + 273.15; %K
h = 10; % W/m^2k


SA = 2 * pi * l_sec_m * d_sec_m; % SA of tube in sun, .5 of normal


Q_sun = 1360 * 0.5 * SA;
%Q_amb = h * SA * (T - T_amb);
Q_rocket = @(T) ((Emissivity * Boltz * T^4 * SA) + (h * SA * (T - T_amb))); 

%Q_test = @(T2) (h * SA * (T_amb - T2))

T_rocket = fzero(@(T) Q_sun-Q_rocket(T),300)


