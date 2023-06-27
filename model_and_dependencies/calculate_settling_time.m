%% set demo defaults
clc; close all; clear all;

set(groot,'defaultAxesFontName','Times')
set(groot,'defaultAxesFontSize',20)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultLineLineWidth', 2);


%% Calculate the settling time

V  = linspace(100,5000,1000); %volume of the melt [km^3] 
d_ice= [50, 80, 100, 150]' ;      %thickness of ice sheet [km]

%Parameters
grav      = 1.352;%gravity: Europa =1.315, Titan = 1.352, [m/s^2]
rho_melt  = 1000; %density of the melt (water) [kg/m^3]
rho_matrix= 917; %density of the bath (ice) [kg/m^3] 
mu        = 1e14;%dynamic viscosity [Pa.s]
yr2s      = 365.25*24*60*60; %year to second conversion [s/year]

%% Calculations 

R    = (3.*V./(4.*pi)).^(1/3);  %calculate radius of the blob [km]

%Calculating the settling velocity

V_settling = (2/9)*(rho_melt - rho_matrix).*(R*1e3).^2.*grav./mu.*(yr2s/1e3); %[km/year]
t_settling = d_ice./ V_settling; %settling time [years]

labels=num2str(d_ice,'$$d_{ice}=$$%d km');

i = 1;
figure(1)
xlim([min(V) max(V)]);
ylabel('Settling time, years');
xlabel('Melt volume, km$$^3$$');
    hold on
while true
    plot(V,t_settling(i,:));
    
    i = i+1;
    if i > length(d_ice), break ; end
end 

legend(labels,'location','best')