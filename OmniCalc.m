%% Constants

R = 8314; % universal gas constant, J/(mol*K)

%Environmental Constants

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol

%Rocket Constants

R_combust = 356; % gas constant of the motor combustion products, J/(kg*K)
T_combust = 2773; % temperature of the motor combustion prodects, K

% I love bees

%% Vehicle Parameters
number_of_bees = 10000;

airframe_diameter = 6; % Airframe Diameter (in)
recovery_bay_length = 30; % Recovery Bay Length (in)
coupler_length = 6; % Coupler Length (in)
shock_cord_length = 400; % Length of Shock Cord (in)

lower_mass = 38; % Lower Section Mass (Ib)
upper_mass = 0.01; % Upper Section Mass (Ib)

drouge_diameter = 36; % Drouge Chute Diameter (in)
packed_drouge_diameter = 2.85; % Packed Drouge Chute Diameter (in)
drouge_n = 2.5; % Drogue Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6)
lower_main_diameter = 144; % Lower Main Chute Diameter (in)
packed_lower_main_diameter = 6; % Packed Lower Main Chute Diameter (in)
lower_main_n = 2.5; % Lower Main Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6)
upper_main_diameter = 144; % Upper Main Chute Diameter (in)
packed_upper_main_diameter = 6; % Packed Upper Main Chute Diameter (in)
upper_main_n = 2.5; % Upper Main Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6)

fin_area = 5.11; % Frontal Fin Area (in^2)

cd_parachute = 2.2; % Parachute Coeffient of Drag
cd_lower = 0.32; % Lower Section Coeffient of Drag
cd_upper = 0.32; % Upper Section Coeffient of Drag
cd_fin = 0.09; % Fin Coeffient of Drag

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
