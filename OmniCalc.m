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
%Ejection: All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment Velocities,

%% Rocket Constants

R_combust = 356; % gas constant of the ejection ejection charge combustion products, J/(kg*K) %%TODO
T_combust = 1833; % temperature of the ejection charge combustion products, K %%TODO

%% Vehicle Parameters

number_of_bees = 10000; % I love bees

upper_mass = 22.5; % Upper Section Mass (Ib)
lower_mass = 44.9; % Lower Section Mass (Ib) %%TODO

% Found from openRocket component analysis
cd_lower = 0.254; % Lower Section Coeffient of Drag (cd_payloadbay + cd_ebay + cd_fincan)/(bottom half of rocket)
cd_upper = 0.109; % Upper Section Coeffient of Drag (cd_upper + cd_nosecone)/(top half of rocket)

internal_volume = 568.09; % Internal volume of the rocket (in^3)

Emissivity = .84; %Assuming a white paint rocket
Length_of_Ebay = 19; %Length of electronics bay in (in) %%Make sure this is the correct section name!!
airframe_outside_diameter= 6.17; %diameter of the rocket (in)

shear_pin_strength = 178; % Tensile strength of shear pins (N) TODO

%% Launch Site Parameters

launch_MSL = 5700; % altitude of the launch site above mean sea level, ft %%TODO
temperature = 91; % ambient temperature of the launch site, F %%TODO
max_wind_vel = 0; % maximum allowable wind speed, (ft/s) %%TODO

%% Constants

R = 8314; % universal gas constant, J/(mol*K)
k_b = 1.38E-23; % Boltzmann Constant (JK^-1)
Boltz = 5.67*10^-8; % Stefan-Boltzmann Constant (W/m^2K^4)
M_a = 0.02897; % Molar Mass of Air (kgmol^-1)
N_A = 6.02E23; % Avagadro's Number (mol^-1)

%% Environmental Constants

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa
rho = 1.225; % atmospheric density at sea level (kg/m^3)

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol

h_amb_air = 10; %heat transfer coeff of ambient air W/m^2*K
ground_wind_speed = 6.5; %wind speed on the ground in m/s

%% Settings

shear_pin_safety_factor = 2; % Safety factor for number of shear pins
is_wind = true; %Is there wind?
dt = 0.01; % Interpolated dt of the Rasaero data (s)
vent_hole_accuracy = 0.0001; % How close internal pressure is to external pressure
vent_hole_presicion = 0.0000254; % How precise the vent holes can be machined (in)

%% Data Imput from RASAero

RASdata = readmatrix('Flight Test2.CSV');
altitudes = (RASdata(:,23).*0.3048)+launch_MSL;
velocities = RASdata(:,18);
velocities_v = RASdata(:,19);
velocities_h = RASdata(:,20);
acclerations = RASdata(:,16);
times = RASdata(:,1);
Ras_dt = RASdata(2,1)-RASdata(1,1);

altitudes_prea = altitudes_to_apogee(altitudes,velocities_v);
velocities_v_posta = velocities_past_apogee(velocities_v);
%altitudes = interpolate_alt(altitudes,Ras_dt,dt);

%% Conversions

upper_mass = upper_mass*4.44822; % Upper Section Mass (N)
lower_mass = lower_mass*4.44822; % Lower Section Mass (N)

launch_MSL = launch_MSL*0.3048; % altitude of the launch site above mean sea level, (m)
temperature = (5/9)*(temperature-32) + 273.15; % ambient temperature of the launch site, (K)
max_wind_vel = max_wind_vel*0.44704; % maximum allowable wind speed, (m/s)

internal_volume = internal_volume*(0.0254^3); % Internal volume of the rocket (m^3)

altitudes_prea = altitudes_prea.*0.3048;
velocities_v_posta = velocities_v_posta.*0.3048;
velocities_h = velocities_h.*0.3048;
acclerations = (acclerations*0.3048)-g;

%% Derived Parameters

total_mass = upper_mass + lower_mass; % Mass of whole rocket (N)
drag_ratio = cd_upper/cd_lower; % The ratio of drag coefficients of the upper and lower airframe

%% Descent Calculations
%% Descent Velocities
[descent_velocity_main,descent_velocity_drogue] = velocity_of_descent(velocities_v_posta);

%% Descent Times
descent_time = time_of_descent(times,altitudes_prea);

%% Landing Kinetic Energies
landing_kinetic_energy = kinetic_energy(velocities(length(velocities)),total_mass);

%% Downrange Drift
downrange_drift = drift(velocities_h,Ras_dt);

%% Ejection Calculations
%All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment 
%% Pre Separation
%% Max Decelartation
max_deceleration = abs(min(acclerations));
max_deceleration = max_deceleration*0.3048;

%% Separation Forces
   % F_upper_lower = drag_force(Cd_fin,rho_max,v_max,A_fin) + drag_force(Cd_lower,rho_max,v_max,A_lower) - drag_force(Cd_upper,rho_max,v_max,A_upper);
F_upper_lower = drag_seperation(max_deceleration,total_mass,drag_ratio,lower_mass);


%% Shear Pins
pins_upper_lower = ceil((F_upper_lower*shear_pin_safety_factor)/shear_pin_strength);

%% Post Separation
%% Ejection Charges
Eject_force = 1.5;

%% Ejection Velocities



%% Parachute Deployment Forces
%% E-Bay Temperature

Ebay_temp = Temp_EBay(ground_wind_speed, Length_of_Ebay, airframe_outside_diameter, Boltz, temperature, h_amb_air, is_wind, Emissivity);

%% Vent Hole

vent_hole_maxTimeSteps = length(altitudes_prea);

t = (0:dt:RASdata(length(altitudes_prea)));

vent_hole_diameter = vent_hole_presicion;
while(true)
    [PRec,P_iRec,P_eRec] = vent_hole_pressure(vent_hole_diameter,dt,vent_hole_maxTimeSteps,P0,internal_volume,k_b,Ebay_temp,temperature,altitudes_prea,rho,launch_MSL,M,g,R/1000,N_A,L);
    if(max(PRec)<vent_hole_accuracy)
        break;
    end
    vent_hole_diameter = vent_hole_diameter + vent_hole_presicion;
end

%% Outputs

fprintf("Descent velocity under main: %3.2fft/s\n",descent_velocity_main/0.3048);
fprintf("Descent velocity under drogue: %3.2fft/s\n",descent_velocity_drogue/0.3048);
fprintf("Descent time: %3.2fs\n",descent_time);
fprintf("Landing kinetic energy: %3.2fNm\n",landing_kinetic_energy);
fprintf("Downrange drift: %3.2fft\n",downrange_drift/0.3048);
fprintf("Maximum decleration: %3.2fft/s^2\n",max_deceleration/0.3048);
fprintf("Number of shear pins for a safety factor of %3.2f: %2.0f\n",shear_pin_safety_factor,pins_upper_lower);
fprintf("Worst case e-bay tempurature on pad: %3.2fF\n",(9/5)*(Ebay_temp-273.15)+32);
fprintf("Minimum vent hole diameter: %1.3fin\n",vent_hole_diameter*39.3701);
fprintf("Internal pressure range: %3.2f-%3.2fpsi\n",max(P_iRec)/6894.76,min(P_iRec(1,length(P_iRec)-1))/6894.76);
fprintf("External pressure range: %3.2f-%3.2fpsi\n",max(P_eRec)/6894.76,min(P_eRec(1,length(P_eRec)-1))/6894.76);

% figure()
% plot(t,P_iRec)
% title("Pressure to Apogee")
% xlabel('Time (s)', 'FontSize', 11)
% ylabel('Internal Pressure (Pa)', 'FontSize', 11)

%% Functions
function altitudes = altitudes_to_apogee(h,v)
    i=1;
    while(v(i)>=0)
        i=i+1;
    end
    altitudes = h(1:i);
end

function velocity = velocities_past_apogee(v)
    i=1;
    while(v(i)>=0)
        i=i+1;
    end
    velocity = v(i:length(v));
end

function altitudes = interpolate_alt(h,dt_current,dt)
    extra_points = round(dt_current/dt);
    altitudes = zeros((length(h)-1)*extra_points,1);
    for(i = 1:length(h)-1)
        for(j = 1:extra_points)
            altitudes(extra_points*(i-1)+j) = h(i) + (j-1)*dt*(h(i+1)-h(i));
        end
    end
end

function [v_descent_main,v_descent_drogue] = velocity_of_descent(velocities)
    v_descent_main = 0;
    v_descent_drogue = 0;
    under_main = false;
    for(i = 2:length(velocities))
        if((velocities(i)-velocities(i-1))/velocities(i)<0.01)
            if(under_main)
                v_descent_main = abs(velocities(i));
            else
                v_descent_drogue = abs(velocities(i));
                under_main = true;
            end
        end
    end
end

function t_descent = time_of_descent(times,altitudes)
    t_a = times(length(altitudes));
    t_f = times(length(times));
    t_descent = t_f - t_a;
end

function KE = kinetic_energy(V,m)
    KE = 0.5*m*V^2;
end

function downrange_drift = drift(velocities_h,dt)
    downrange_drift = sum(velocities_h.*dt);
end

function Fsep = drag_seperation(a,M,R,M1)
    Fsep = a*(M/(1+R)-M1);
end

function D = drag_force(Cd,rho,v,A)
    D = (1/2)*Cd*rho*v^2*A;
end

function T_EBay = Temp_EBay(ground_wind_speed, Length_of_Ebay, airframe_outside_diameter, Boltz, temperature, h_amb_air, is_wind, Emissivity) %Worst case temp of EBay, assuming no wind

    if is_wind
        h_forced = 10.45 - ground_wind_speed + 10*sqrt(ground_wind_speed); %This is for the cooling from wind, I've found different sources giving different coefficents,
        SA_Ebay = 2 * pi * Length_of_Ebay * airframe_outside_diameter; % Surface Area of Ebay
        Q_sun = 1360 * 0.5 * SA_Ebay;
        Q_EBay = @(T) ((Emissivity * Boltz * T^4 * SA_Ebay) + (h_forced * 0.5*SA_Ebay * (T - temperature))); %Assuming only half of the Surface recieves wind
        T_EBay = fzero(@(T) Q_sun-Q_EBay(T),300);
    else
        SA_Ebay = 2 * pi * Length_of_Ebay * airframe_outside_diameter; % Surface Area of Ebay
        Q_sun = 1360 * 0.5 * SA_Ebay;
        Q_EBay = @(T) ((Emissivity * Boltz * T^4 * SA_Ebay) + (h_amb_air * SA_Ebay * (T - temperature))); %for worst case scenario assume no wind cooling
        T_EBay = fzero(@(T) Q_sun-Q_EBay(T),300);
    end

end 

function [PRec,P_iRec,P_eRec] = vent_hole_pressure(d,dt,maxTimeSteps,P_0,V,k_b,T_0,Tout_0,altitudes,rho,h_0,M_a,g,R,N_A,L)

    A = pi*(d/2)^2;
    P_out = pressure_at_altitude(P_0,g,M_a,altitudes(1),0,R,Tout_0);
    N = (P_out*V)/(k_b*T_0);
    xCurr = [N,altitudes(1)]; 
    PRec = zeros(1,maxTimeSteps);
    PRec(1) = 0;

    P_iRec = zeros(1,maxTimeSteps);
    P_iRec(1) = P_out;

    P_eRec = zeros(1,maxTimeSteps);
    P_eRec(1) = P_out;
    
    i=2;
    while(i<length(altitudes))

        k1=vent_hole_dynamics(xCurr,V,M_a,rho,T_0,Tout_0,P_0,k_b,A,h_0,R,g,N_A,L,dt,altitudes(i),altitudes(i+1))*dt;
        k2=vent_hole_dynamics(xCurr+1/2*k1,V,M_a,rho,T_0,Tout_0,P_0,k_b,A,h_0,R,g,N_A,L,dt,altitudes(i),altitudes(i+1))*dt;
        k3=vent_hole_dynamics(xCurr+1/2*k2,V,M_a,rho,T_0,Tout_0,P_0,k_b,A,h_0,R,g,N_A,L,dt,altitudes(i),altitudes(i+1))*dt;
        k4=vent_hole_dynamics(xCurr+k3,V,M_a,rho,T_0,Tout_0,P_0,k_b,A,h_0,R,g,N_A,L,dt,altitudes(i),altitudes(i+1))*dt;
        xCurr=xCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;
        
        T = tempurature_at_altitude(Tout_0,altitudes(i),h_0,L);
        P_out = pressure_at_altitude(P_0,g,M_a,altitudes(i+1),0,R,T);

        N = xCurr(1);
        P_in = (N/V)*k_b*T_0;
    
        PRec(i) = (P_in-P_out)/P_out;
        P_iRec(i) = P_in;
        P_eRec(i) = P_out;
        i=i+1;
    end
    P_iRec(length(P_iRec)) = P_iRec(length(P_iRec)-1);
    P_eRec(length(P_eRec)) = P_eRec(length(P_eRec)-1);
end


function xDot = vent_hole_dynamics(x,V,M_a,rho,T_0,Tout_0,P_0,k_b,A,h_0,R,g,N_A,L,dt,altCurr,altNext)
    N = x(1); h = x(2);
    T = tempurature_at_altitude(Tout_0,h,h_0,L);
    P_out = pressure_at_altitude(P_0,g,M_a,h,0,R,T);

    P_in = (N/V)*k_b*T_0;
    mDot = mass_flow_rate(A,P_in,P_out,rho);

    NDot = -(N_A/M_a)*mDot;
    hDot = (altNext-altCurr)/dt;
    xDot = [NDot,hDot];
end

function mDot = mass_flow_rate(A,P_in,P_out,rho)
        mDot = ((P_in-P_out)/abs(P_in-P_out))*A*sqrt(2*rho*abs(P_in-P_out));
end

function T = tempurature_at_altitude(T_0,h,h_0,L)
    T = T_0 - L*(h-h_0);
end

function P = pressure_at_altitude(P_0,g,M,h,h_0,R,T)
    P = P_0*exp((-g*M*(h-h_0))/(R*T));
end

