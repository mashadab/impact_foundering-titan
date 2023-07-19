%Analytical 
set(groot, 'DefaultFigureVisible', 'on');
clc; close all; clear all;
%% Input variables
%Thermodynamics
Ti =  223.16;       %surrounding ice temperature [K]

%Settling
%{
%Titan
grav      = 1.352;%gravity: Europa =1.315, Titan = 1.352, [m/s^2]
V  = linspace(100,5000,1000); %volume of the melt [km^3] 
d_ice= [50, 80, 100, 150]' ;  %thickness of ice sheet [km]
%}

%Europa
grav= 1.315; %Titan = 1.352, [m/s^2]
V  = 31%linspace(0.1,1000,1000); %volume of the melt [km^3] 
d_ice= [10]' ;  %thickness of ice sheet [km]
% Calculations 
R    = (3.*V./(4.*pi)).^(1/3);  %calculate radius of the blob [km]
yr2s = 365.25*24*60*60; %year to second conversion [s/year]

%% Thermodynamic parameters
L    = 333.6e3;  %latent heat of fusion [J/kg]
cp_i = 2.1e3;    %Specific heat of ice [J/kg-K]
rho_i= 917;      %density of ice [kg/m^3]
Tm =  273.16;    %melting temperature; temperature of the drop [K]
k  =  2.22;       %thermal conductivity [W/m-K]
kappa = k / (rho_i * cp_i); %Thermal diffusivity [m^2/s]



%% Calculate solidification: Stefan problem
tau_max = @(beta)  1./6 + 1./(6.*beta) - 1./(3.*(2.*pi).^0.5.*beta.^(1.5));


figure();
beta_array = logspace(-1,3,1000);
loglog(beta_array,tau_max(beta_array))
xlabel 'Beta'
ylabel 'Dim-less time, Tau'


beta = L/(cp_i*(Tm - Ti)); %Stefan number

S = linspace(0,1-1e-10,10000); %S is the dimless location of the front w.r.t initial position
tau  =  (3.*S.^2 - 2.*S.^3)./6 + S.^2 ./ (6.*beta) - S.^2./(45.*beta.^2.*(1-S)); %dimless time

figure();
plot(tau,S);
xlim([0 0.3])
ylabel 'S'
xlabel 'Dimless time, tau'
t_solid = tau_max(beta) * beta * (R*1e3).^2 / (kappa * yr2s) %redimensionalizing time

%% Calculate Settling: Stokes settling
%Parameters
rho_melt  = 1000; %density of the melt (water) [kg/m^3]
rho_matrix= 917;  %density of the bath (ice) [kg/m^3] 
mu        = 1e14; %dynamic viscosity [Pa.s]
yr2s      = 365.25*24*60*60; %year to second conversion [s/year]

%Calculating the settling velocity

V_settling = (2/9)*(rho_melt - rho_matrix).*(R*1e3).^2.*grav./mu.*(yr2s/1e3); %[km/year]
t_settling_without_refreezing = d_ice./ V_settling %settling time [years]
blacks = linspace(0.9,0,length(d_ice));

labels=num2str(d_ice,'$$d_{ice}=$$%d km');


%Variable radius settling

Var_R = (1 - S).*R   %Transient radius of the drop
V_settling_var = (2/9)*(rho_melt - rho_matrix).*(Var_R*1e3).^2.*grav./mu.*(yr2s/1e3); %[km/year]
timestamp_var  = tau .* beta * (R*1e3).^2 ./ (kappa * yr2s) %redimensionalizing time array [years]

figure();
loglog(timestamp_var,V_settling_var)
ylabel 'V [km/year]'
xlabel 'Time [year]'

i = 1;
z = zeros(1,9999);
while true
    
    depth_var(i) = V_settling_var(i).*(timestamp_var(i+1) - timestamp_var(i)); %instantaneous depth [km]
    
    i = i+1
    
    if i >= 10000, break ; end
    
end 
    

figure();
loglog(timestamp_var(:,1:end-1),depth_var)
xlabel 'T [years]'
ylabel 'Depth [km]'