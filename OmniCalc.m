clear; close all; clc;

%To do from George: 
%Double check units on everything, all inputs and outputs in imperial 
%Test with same CFD numbers (500ft/s, plus whatever pressure, density, etc. changes with alt) 
%General To do: 
%Improve the program: 
%Use functions everywhere possible
%Remove anything trajectory related 
%Re-check equations: 
%Re-derive equations used 
%Other 
%Keep velocities at all stages of the parachute 

%Livescript in conversion factors - do this ish later tho - T Champ
%Descent Velocity: Get certain values from Rasaero data
%Descent Time: Rasaero
%Landuing Kinetic Energy: We calc
%Downrange drft: Rasaero
%Ejection: All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment Velocities,
%Cut all the function stuff at the bottom

%% Data Imput from RASAero

data_table_acc = readtable('.csv');
%% Rocket Constants

R_combust = 356; % gas constant of the motor combustion products, J/(kg*K) %%TODO
T_combust = 1833; % temperature of the motor combustion prodects, K %%TODO

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

drag_ratio = 0; % the ratio of drag coefficients of the upper and lower airframe %%TODO

internal_volume = 568.09; % Internal volume of the rocket (in^3)

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
k_b = 1.38E-23; % Boltzman Constant (JK^-1)
M_a = 0.02897; % Molar Mass of Air (kgmol^-1)


%% Environmental Constants

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol


%% Settings

vent_hole_accuracy = 0.001; % How close internal pressure is to external pressure
vent_hole_presicion = 0.01; % How precise the vent holes can be machined 
vent_hole_dt = 0.01; % dt for vent hole calculation (s)
vent_hole_maxTimeSteps = 300; % Max time for vent hole calculation (s)


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

internal_volume = internal_volume*(0.0254^3); % Internal volume of the rocket (m^3)


%% Rasaero Integration




%% Descent Calculations
%% Descent Velocities
%Get certain values from Rasaero data

%% Descent Times
%Rasaero

%% Landing Kinetic Energies
%We calc

%% Downrange Drift
%Rasaero

%% Ejection Calculations
%All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment 
%% Pre Separation
%% Separation Forces
   % F_upper_lower = drag_force(Cd_fin,rho_max,v_max,A_fin) + drag_force(Cd_lower,rho_max,v_max,A_lower) - drag_force(Cd_upper,rho_max,v_max,A_upper);
    % F_lower_fin = drag_force(Cd_fin,rho_max,v_max,A_fin) - drag_force(Cd_lower,rho_max,v_max,A_lower) - drag_force(Cd_upper,rho_max,v_max,A_upper);

%% Shear Pins
    pins_upper_lower = ceil((F_upper_lower*1.25)/shear_pin_strength);
    pins_lower_fins = ceil((F_upper_lower*1.25)/shear_pin_strength);

%% Post Separation
%% Ejection Charges
    Eject_force = 1.5*

%% Ejection Velocities


%% Parachute Deployment Forces


%% Vent Hole

vent_hole_diameter = vent_hole_presicion; % First Guess for vent hole diameter
while(true)
    PRec = vent_hole_pressure(vent_hole_dt,vent_hole_maxTimeSteps,P0,internal_volume,k_b,internal_temperature,altitudes,launch_MSL,M);
    if(max(PRec)<vent_hole_accuracy)
        break;
    end
    vent_hole_diameter = vent_hole_diameter + vent_hole_presicion;
end

%% Functions
function D = drag_force(Cd,rho,v,A)
    D = (1/2)*Cd*rho*v^2*A;
    
end

function PRec = vent_hole_pressure(dt,maxTime,P_0,V,k_b,T_0,altitude,h_0,rho,M_a)

    t = zeros(1,maxTime/dt);
    
    N = (P_0*V)/(k_b*T_0);
    xCurr = N; 
    PRec = zeros(1,length(t));
    PRec(1) = P_0;
    
    i=2;
    while(true)
        
        T = tempurature_at_altitude(T_0,altitude(i),h_0,rho);
        P_out = pressure_at_altitude(P_0,g,M,altitude(i),h_0,R,T);
        k1=vent_hole_dynamics(xCurr,V,P_out,M_a,rho,T_0)*dt;
        k2=vent_hole_dynamics(xCurr+1/2*k1,V,P_out,M_a,rho,T_0)*dt;
        k3=vent_hole_dynamics(xCurr+1/2*k2,V,P_out,M_a,rho,T_0)*dt;
        k4=vent_hole_dynamics(xCurr+k3,V,P_out,M_a,rho,T_0)*dt;
        xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
        
        t(i)=t(i-1)+dt;

        N = xCurr;
        P = (N/V)*k_b*T_0;
    
        PRec(i) = abs((P-P_0)/P_0);
        i=i+1;
    end
end


function xDot = vent_hole_dynamics(N,V,P_out,M_a,rho,T_0)
    P_in = (N/V)*k_b*T_0;
    mDot = mass_flow_rate(A,P_in,P_out,rho);

    NDot = -(N/M_a)*mDot;
    xDot = NDot;
end

function mDot = mass_flow_rate(A,P_in,P_out,rho)
    mDot = A*sqrt(2*rho*(P_in-P_out));
end

function T = tempurature_at_altitude(T_0,h,h_0,L)
    T = T_0 + L*(h-h_0);
end

function P = pressure_at_altitude(P_0,g,M,h,h_0,R,T)
    P = P_0*exp((-g*M*(h-h_0))/(R*T));
end

function KE = Kinetic_E(V,m)
KE = 0.5*m*V^2;
end

