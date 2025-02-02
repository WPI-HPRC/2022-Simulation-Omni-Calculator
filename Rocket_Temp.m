%% Comments
% Possibly the h's should be combined? They refer to stagnant and moving
% air, but in reality it's impossible to have both at once
%Should we assume that wind can only hit one side of the rocket at a time?
%If we're looking for a worst case scenario we should assume zero wind.


Boltz = 5.67*10^-8; % W/m^2K^4
Emissivity = .84;
l_sec = 19; %in
l_sec_m = l_sec * 0.0254; %m
d_sec = 6.17; %in
d_sec_m = d_sec * 0.0254; %m
T_amb_F = 91; %F
T_amb = (T_amb_F - 32) * 5/9 + 273.15; %K
h = 10; % W/m^2k ambient air
v = 6.5; %Wind speed in m/s

%h_forced = 10.45 - v + 10*sqrt(v); %This is for the cooling from wind, I've found different sources giving different coefficents

SA = 2 * pi * l_sec_m * d_sec_m; % SA of tube in sun, .5 of normal


Q_sun = 1360 * 0.5 * SA;
%Q_amb = h * SA * (T - T_amb); %This is really just the cooling if there
%was no wind
Q_rocket = @(T) ((Emissivity * Boltz * T^4 * SA) + (h * SA * (T - T_amb)));  

%Q_wind = h_forced * .5*SA *(T-T_amb)

T_rocket = fzero(@(T) Q_sun-Q_rocket(T),300);

T_Rocket_F = ( T_rocket - 273.15) * 9/5 + 32


