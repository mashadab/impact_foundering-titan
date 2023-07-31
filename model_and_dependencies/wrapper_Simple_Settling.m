%% demo for manannan crater simulation in figure 3
clear all; close all;
%% model parameters
cox_and_bauer_impact_key = 'cold1';
basal_viscosity = 1e14; % basal viscosity of ice
viscosity_activation_energy = 50e3; % viscosity activation energy
%% model run simulation
driver_func_Simple_Settling(cox_and_bauer_impact_key,basal_viscosity, ...
    viscosity_activation_energy) %main driver routine
         