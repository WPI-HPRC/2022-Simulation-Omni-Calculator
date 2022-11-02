%% Constants

R = 8314; % universal gas constant, J/(mol*K)

%Environmental Constants

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol

%Rocket Constants

R_combust = 119.2; % gas constant of the motor combustion products, J/(kg*K)
T_combust = 1837.2; % temperature of the motor combustion prodects, K

% I love bees

%% Vehicle Parameters
number_of_bees = 10000;

%% Flight parameters
burnout_AGL = 1477; % predicted burnout altitude above ground level, ft
apogee_AGL = 4099; % predicted highest point of flight above ground level, ft
main_AGL = 600; % predectied altitude above ground level, ft

Max_Vel = 58; % the predicted maximum velocity of the rocket, ft
Max_drift = 2500; % maximum allowable drift, ft

%% Launch Site Parameters
launch_MSL = 900; % altitude of the launch site above mean sea level, ft
temperature = 70; % ambient temperature of the launch site, F
max_wind_vel = 20; % maximum allowable wind speed, mph
