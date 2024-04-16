%% demo for manannan crater simulation in figure 3
clear all; close all;
%% model parameters
% manannan crater key from Cox and Bauer, 2015 simulations
cox_and_bauer_impact_key = 'Wakita-Stokes';
basal_viscosity = 1e14; % basal viscosity of ice
viscosity_activation_energy = 50e3; % viscosity activation energy
%% model run simulation
impactorTempMeltFuncDSModPresForm2StokesTitandataShigeruTest(cox_and_bauer_impact_key,basal_viscosity, ...
    viscosity_activation_energy) %main driver routine

%impactorTempMeltFuncTitan_data : original Stokes code with Titan data
%impactorTempMeltFuncDarcyStokesModPresForm2StokesTitandata : new
%impactorTempMeltFuncDSModPresForm2StokesTitandataShigeru: new Shigeru data
%impactorTempMeltFuncDSModPresForm2StokesTitandataShigeruTest: new Shigeru data for testing
%simpleporositywavetest: Porosity wave test
%Darcy-Stokes with Titan iSALE data and kc handle
         