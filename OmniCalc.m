
%% Rocket Constants

R_combust = 356; % gas constant of the motor combustion products, J/(kg*K) %%TODO
T_combust = 2773; % temperature of the motor combustion prodects, K %%TODO

% I love bees


%% Vehicle Parameters

number_of_bees = 10000;

airframe_diameter = 6; % Airframe Diameter (in) %%TODO
recovery_bay_length = 30; % Recovery Bay Length (in) %%TODO
coupler_length = 6; % Coupler Length (in) %%TODO
shock_cord_length = 400; % Length of Shock Cord (in) %%TODO

lower_mass = 38; % Lower Section Mass (Ib) %%TODO
upper_mass = 0.01; % Upper Section Mass (Ib) %%TODO

drouge_diameter = 36; % Drouge Chute Diameter (in) %%TODO
packed_drouge_diameter = 2.85; % Packed Drouge Chute Diameter (in) %%TODO
drouge_n = 2.5; % Drogue Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6) %%TODO
lower_main_diameter = 144; % Lower Main Chute Diameter (in) %%TODO
packed_lower_main_diameter = 6; % Packed Lower Main Chute Diameter (in) %%TODO
lower_main_n = 2.5; % Lower Main Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6) %%TODO
upper_main_diameter = 144; % Upper Main Chute Diameter (in) %%TODO
packed_upper_main_diameter = 6; % Packed Upper Main Chute Diameter (in) %%TODO
upper_main_n = 2.5; % Upper Main Chute Canopy Fill Constant (see Parachute Recovery Systems Design Manual, Knacke, Table 5-6) %%TODO

fin_area = 5.11; % Frontal Fin Area (in^2) %%TODO

cd_parachute = 2.2; % Parachute Coeffient of Drag %%TODO
cd_lower = 0.32; % Lower Section Coeffient of Drag %%TODO
cd_upper = 0.32; % Upper Section Coeffient of Drag %%TODO
cd_fin = 0.09; % Fin Coeffient of Drag %%TODO


%% Flight parameters

burnout_AGL = 1477; % predicted burnout altitude above ground level, ft %%TODO
apogee_AGL = 4099; % predicted highest point of flight above ground level, ft %%TODO
main_AGL = 600; % predectied altitude above ground level, ft %%TODO

Max_Vel = 58; % the predicted maximum velocity of the rocket, ft %%TODO
Max_drift = 2500; % maximum allowable drift, ft %%TODO


%% Launch Site Parameters

launch_MSL = 5700; % altitude of the launch site above mean sea level, ft %%TODO
temperature = 90; % ambient temperature of the launch site, F %%TODO
max_wind_vel = 20; % maximum allowable wind speed, mph %%TODO


%% Constants

R = 8314; % universal gas constant, J/(mol*K)


%% Environmental Constants

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol


%% Conversions

airframe_diameter = airframe_diameter*0.0254; % Airframe Diameter (m)
recovery_bay_length = recovery_bay_length*0.0254; % Recovery Bay Length (m)
coupler_length = coupler_length*0.0254; % Coupler Length (m)
shock_cord_length = shock_cord_length*0.0254; % Length of Shock Cord (m)

lower_mass = lower_mass*4.44822; % Lower Section Mass (N)
upper_mass = upper_mass*4.44822; % Upper Section Mass (N)

drouge_diameter = drouge_diameter*0.0254; % Drouge Chute Diameter (m)
packed_drouge_diameter = packed_drouge_diameter*0.0254; % Packed Drouge Chute Diameter (m)
lower_main_diameter = lower_main_diameter*0.0254; % Lower Main Chute Diameter (m)
packed_lower_main_diameter = packed_lower_main_diameter*0.0254; % Packed Lower Main Chute Diameter (m)
upper_main_diameter = upper_main_diameter*0.0254; % Upper Main Chute Diameter (m)
packed_upper_main_diameter = packed_upper_main_diameter*0.0254; % Packed Upper Main Chute Diameter (m)

fin_area = fin_area*0.0254*0.0254; % Frontal Fin Area (m^2)

burnout_AGL = burnout_AGL*0.3048; % predicted burnout altitude above ground level, (m)
apogee_AGL = apogee_AGL*0.3048; % predicted highest point of flight above ground level, (m)
main_AGL = main_AGL*0.3048; % predectied altitude above ground level, (m)

Max_Vel = Max_Vel*0.3048; % the predicted maximum velocity of the rocket, (m)
Max_drift = Max_drift*0.3048; % maximum allowable drift, (m)

launch_MSL = launch_MSL*0.3048; % altitude of the launch site above mean sea level, (m)
temperature = (5/9)*(temperature-32); % ambient temperature of the launch site, (C)
max_wind_vel = max_wind_vel*0.44704; % maximum allowable wind speed, (m/s)



