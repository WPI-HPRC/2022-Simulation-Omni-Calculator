clear; close all; clc;

%To do from George: 
%Double check units on everything, all inputs and outputs in imperial 
%Test with same CFD numbers (500ft/s, plus whatever pressure, density, etc. changes with alt) 
%General To do: 
%Cl vs alpha graph (verify finsim results)
%Base drag (if we run out of things to do)
%Livescript in conversion factors - do this ish later tho - T Champ
%Ejection: All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment Velocities,

%% Rocket Constants

R_combust = 356; % gas constant of the ejection ejection charge combustion products, J/(kg*K) %%TODO
T_combust = 1833; % temperature of the ejection charge combustion products, K %%TODO

%% Vehicle Parameters

number_of_bees = 10000; % I love bees

upper_mass = 22.5; % Upper Section Mass (Ib)
lower_mass = 44.9; % Lower Section Mass (Ib) %%TODO
Rocket_Surface_Area = 4*(91.11 + 531.65 + 26.02); % Total surfaec area of Rocket body (no fins)

% Found from openRocket component analysis
cd_lower = 0.254; % Lower Section Coeffient of Drag (cd_payloadbay + cd_ebay + cd_fincan)/(bottom half of rocket)
cd_upper = 0.109; % Upper Section Coeffient of Drag (cd_upper + cd_nosecone)/(top half of rocket)

internal_volume = 568.09; % Internal volume of the rocket (in^3)

Emissivity = .84; %Assuming a white paint rocket
Length_of_Ebay = 19; %Length of electronics bay in (in) %%Make sure this is the correct section name!!
airframe_outside_diameter= 6.17; %diameter of the rocket (in)

shear_pin_strength = 178; % Tensile strength of shear pins (N) TODO

%% Constants

R = 8314; % universal gas constant, J/(mol*K)
k_b = 1.38E-23; % Boltzmann Constant (JK^-1)
sigma = 5.67*10^-8; % Stefan-Boltzmann Constant (W/m^2K^4)
M_a = 0.02897; % Molar Mass of Air (kgmol^-1)
N_A = 6.02E23; % Avagadro's Number (mol^-1)

%% Environmental Constants

launch_MSL = 5700; % altitude of the launch site above mean sea level, ft %%TODO
temperature = 110; % ambient temperature of the launch site, F %%TODO

g = 9.81; % acceleration due to Earth's gravity, m/s    
P0 = 101325; % atmospheric pressure at sea level, Pa
rho = 1.225; % atmospheric density at sea level (kg/m^3)
vol = 42069; % volume of air (m^3)

R_air = 287.058; % gas constant of air, J/(kg*K)
L = 0.00976; % temperature Lapse rate of air, K/m
M = 0.02896968; % molar mass of air, kg/mol
Ca_air = .42069; %  SpHeat Air
M_air = .42069; %kg   Mass air

h_amb_air = 10; %heat transfer coeff of ambient air W/m^2*K
ground_wind_speed = 6.5; %wind speed on the ground in m/s

Sun = 1360; %w/m^2   Heat from the sun

%% Settings

shear_pin_safety_factor = 2; % Safety factor for number of shear pins
is_wind = false; %Is there wind?
dt = 0.01; % Interpolated dt of the Rasaero data (s)
vent_hole_accuracy = 0.001; % How close internal pressure is to external pressure
vent_hole_presicion = 0.0000254; % How precise the vent holes can be machined (in)

%% Data Imput from RASAero

RASdata = readmatrix('Flight Test2.CSV');
altitudes = (RASdata(:,23).*0.3048)+launch_MSL*0.3048;
velocities = RASdata(:,18);
velocities_v = RASdata(:,19);
velocities_h = RASdata(:,20);
acclerations = RASdata(:,16);
times = RASdata(:,1);
Ras_dt = RASdata(2,1)-RASdata(1,1);

altitudes_prea = varriable_to_apogee(altitudes,velocities_v);
velocities_v_posta = varriable_past_apogee(velocities_v,velocities_v);
acclerations_posta = varriable_past_apogee(acclerations,velocities_v);
%altitudes = interpolate_alt(altitudes,Ras_dt,dt);

%% Conversions

upper_mass = upper_mass*4.44822; % Upper Section Mass (N)
lower_mass = lower_mass*4.44822; % Lower Section Mass (N)

Length_of_Ebay = Length_of_Ebay*0.0254;
airframe_outside_diameter = airframe_outside_diameter*0.0254;

launch_MSL = launch_MSL*0.3048; % altitude of the launch site above mean sea level, (m)
temperature = (5/9)*(temperature-32) + 273.15; % ambient temperature of the launch site, (K)

internal_volume = internal_volume*(0.0254^3); % Internal volume of the rocket (m^3)

Rocket_Surface_Area = Rocket_Surface_Area*0.00064516;

altitudes_prea = altitudes_prea.*0.3048;
velocities_v_posta = velocities_v_posta.*0.3048;
velocities_h = velocities_h.*0.3048;
acclerations = (acclerations*0.3048);

%% Derived Parameters

total_mass = upper_mass + lower_mass; % Mass of whole rocket (N)
drag_ratio = cd_upper/cd_lower; % The ratio of drag coefficients of the upper and lower airframe

%% Descent Calculations
%% Descent Velocities
[descent_velocity_main,descent_velocity_drogue] = velocity_of_descent(velocities_v_posta,acclerations_posta);

%% Descent Times
descent_time = time_of_descent(times,altitudes_prea);

%% Landing Kinetic Energies
upper_landing_kinetic_energy = kinetic_energy(velocities(length(velocities)),upper_mass);
lower_landing_kinetic_energy = kinetic_energy(velocities(length(velocities)),lower_mass);

%% Downrange Drift
downrange_drift = drift(velocities_h,Ras_dt);

%% Ejection Calculations
%All us including, Separation Forces, Shear Pins, Ejection Charges, Ejection Velocities, Parachute Deployment 
%% Pre Separation
%% Max Decelartation
max_deceleration = abs(min(acclerations+g));

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

Ebay_temp = Temp_EBay(ground_wind_speed, Length_of_Ebay, airframe_outside_diameter, sigma, temperature, h_amb_air, is_wind, Emissivity,Sun);

%% Rocket Temp over time

TCurr = temperature;
TRec = zeros(1,length(velocities));
TRec(1) = TCurr;

i = 2;
while(i < length(velocities))
    T_o = tempurature_at_altitude(temperature,altitudes(i),launch_MSL,L);
    k1=temperature_dynamics(TCurr,T_o,ground_wind_speed,velocities(i),Rocket_Surface_Area,Emissivity,sigma,Sun,Ca_air,M_air)*dt;
    k2=temperature_dynamics(TCurr+1/2*k1,T_o,ground_wind_speed,velocities(i),Rocket_Surface_Area,Emissivity,sigma,Sun,Ca_air,M_air)*dt;
    k3=temperature_dynamics(TCurr+1/2*k2,T_o,ground_wind_speed,velocities(i),Rocket_Surface_Area,Emissivity,sigma,Sun,Ca_air,M_air)*dt;
    k4=temperature_dynamics(TCurr+k3,T_o,ground_wind_speed,velocities(i),Rocket_Surface_Area,Emissivity,sigma,Sun,Ca_air,M_air)*dt;
    TCurr=TCurr+1/6*k1+1/3*k2+1/3*k3+1/6*k4;

    TRec(i) = TCurr;
    i = i+1;
end

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



% Density

density_prea = zeros(length(altitudes_prea),1);
for(i = 1:length(altitudes_prea))
    density_prea(i) = density_at_altitude(P0,g,M,altitudes_prea(i),0,R/1000,temperature,L);
end

%% Outputs

fprintf("Descent velocity under main: %3.2fft/s\n",descent_velocity_main/0.3048);
fprintf("Descent velocity under drogue: %3.2fft/s\n",descent_velocity_drogue/0.3048);

fprintf("Descent time: %3.2fs\n",descent_time);

fprintf("Landing kinetic energy of upper section: %3.2fNm\n",upper_landing_kinetic_energy);
fprintf("Landing kinetic energy of lower section: %3.2fNm\n",lower_landing_kinetic_energy);

fprintf("Downrange drift: %3.2fft\n",downrange_drift/0.3048);

fprintf("Maximum decleration: %3.2fft/s^2\n",max_deceleration/0.3048);

fprintf("Speration Force: %3.2fN\n",F_upper_lower);
fprintf("Number of shear pins for a safety factor of %3.2f: %2.0f\n",shear_pin_safety_factor,pins_upper_lower);
fprintf("Worst case e-bay tempurature on pad: %3.2fF\n",(9/5)*(Ebay_temp-273.15)+32);
fprintf("Min e-bay tempurature during flight: %3.2fF\n",((9/5)*((min(TRec(1:length(TRec)-1)))-273.15)+32));

fprintf("Minimum vent hole diameter: %1.3fin\n",vent_hole_diameter*39.3701);
fprintf("Internal pressure range: %3.2f-%3.2fpsi\n",max(P_iRec)/6894.76,min(P_iRec(1,length(P_iRec)-1))/6894.76);
fprintf("External pressure range: %3.2f-%3.2fpsi\n",max(P_eRec)/6894.76,min(P_eRec(1,length(P_eRec)-1))/6894.76);

% figure()
% plot(t,P_iRec)
% title("Pressure to Apogee")
% xlabel('Time (s)', 'FontSize', 11)
% ylabel('Internal Pressure (Pa)', 'FontSize', 11)

figure(2)
%plot(times(1:1:length(TRec)-1),(9/5)*(TRec(1:length(TRec)-1)-273.15)+32)
title("Temp over Flight")
xlabel('Time (s)', 'FontSize', 11)
ylabel('Temp (K)', 'FontSize', 11)


Ck = (0.981*((7.3*2.2)^(3/2)))/29
Fsep_drogue = ((.5*.981*108^2)*(0.66*0.97)*1)/21
Fsep_main = ((.5*.9877*28.4^2)*(7.3*2.2)*1)/21




%% Functions
function new_varriable = varriable_to_apogee(varriable,velocities)
    i=1;
    while(velocities(i)>=0)
        i=i+1;
    end
    new_varriable = varriable(1:i);
end

function new_varriable = varriable_past_apogee(varriable,velocities)
    i=1;
    while(velocities(i)>=0)
        i=i+1;
    end
    new_varriable = varriable(i:length(velocities));
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

function [v_descent_main,v_descent_drogue] = velocity_of_descent(v,a)
    largest_a_jump = 0;
    v_descent_drogue = 0;
    for(i = 2:length(v))
        if(abs(a(i)-a(i-1))>largest_a_jump)
            largest_a_jump = abs(a(i)-a(i-1));
            v_descent_drogue = abs(v(i-1));
        end
    end
    v_descent_main = abs(v(length(v)));
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

function T_EBay = Temp_EBay(ground_wind_speed, Length_of_Ebay, airframe_outside_diameter, Boltz, temperature, h_amb_air, is_wind, Emissivity,Sun) %Worst case temp of EBay, assuming no wind

    if is_wind
        h_forced = 12.12 - 1.16*ground_wind_speed + 11.6*sqrt(ground_wind_speed); %This is for the cooling from wind, I've found different sources giving different coefficents,
        SA_Ebay = 2 * pi * Length_of_Ebay * airframe_outside_diameter; % Surface Area of Ebay
        Q_sun = Sun * 0.5 * SA_Ebay;
        Q_EBay = @(T) ((Emissivity * Boltz * (T^4) * SA_Ebay) + (h_forced * 0.5*SA_Ebay * (T - temperature))); %Assuming only half of the Surface recieves wind
        T_EBay = fzero(@(T) Q_sun-Q_EBay(T),300);
    else
        SA_Ebay = 2 * pi * Length_of_Ebay * airframe_outside_diameter; % Surface Area of Ebay
        Q_sun = 1360 * 0.5 * SA_Ebay;
        Q_EBay = @(T) ((Emissivity * Boltz * (T^4) * SA_Ebay) + (h_amb_air * SA_Ebay * (T - temperature))); %for worst case scenario assume no wind cooling
        T_EBay = fzero(@(T) Q_sun-Q_EBay(T),300);
    end

end 

function T_dot = temperature_dynamics(T,temperature,ground_wind_speed,velocitiy,Rocket_Surface_Area,Emissivity,sigma,Sun,Ca_air,M_air)
    % Y is up, X is horizontal
    h_forced_x = 12.12 - 1.16*ground_wind_speed + 11.6*sqrt(ground_wind_speed); %This is for the cooling from wind
    h_forced_y = 12.12 - 1.16*velocitiy + 11.6*sqrt(velocitiy); %This is for the cooling from wind
    Half_Rocket_SA = Rocket_Surface_Area/2; % Surface Area of Ebay
    Area_Rocket_Y = pi * (.0254/2)^2; %m^2   Z axis Area of the rocket 
    Q_sun = Sun * Area_Rocket_Y;%   Heat from the sun applied to rocket
    Q_Y_cool = h_forced_y * Area_Rocket_Y * (T - temperature);
    Q_X_cool = h_forced_x * Half_Rocket_SA * (T - temperature);
    Qrad = sigma * Emissivity *Rocket_Surface_Area * (T^4 - temperature^4);
    Q = Q_sun - (Q_Y_cool + Q_X_cool + Qrad);
    T_dot = Q/(Ca_air * M_air);
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

function rho = density_at_altitude(P_0,g,M,h,h_0,R,T_0,L) %%TODO
    T = tempurature_at_altitude(T_0,h,h_0,L);
    P = pressure_at_altitude(P_0,g,M,h,h_0,R,T);
    rho = P/(R*T);
end

