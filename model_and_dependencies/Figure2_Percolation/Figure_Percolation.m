% file: Figure1_Overview.m
% author: Marc Hesse
% description: Introductory figure for oxidant percolation paper.
% date: 21 May 2020

clear all,close all, clc
% set(groot,'ScreenPixelsPerInch',218); 
set(0, 'DefaultFigureRenderer', 'painters'); 
col = marc_colors();

%set(groot,'defaultAxesFontName','Times')
%set(groot,'defaultAxesFontSize',20)
%set(groot,'defaulttextinterpreter','latex')
%set(groot,'defaultAxesTickLabelInterpreter','latex')
%set(groot,'defaultLegendInterpreter','latex')
%set(groot, 'DefaultFigureVisible', 'on');

% Figue size
scsz = get(0,'ScreenSize');
scw = scsz(3); sch = scsz(4);

specs.linewidth = 1.0;
specs.fontsize_axis = 18;
specs.fontsize_axis_label = 18;
specs.fontsize_label = 18;
specs.fontsize_text = 18;
specs.fontsize_subpanel = 18;

% Nature figure pixels
% single column: 1040 pixels
% double column: 2080 pixels
% figW = 600/2; %scw/6;
figW = 18;

gapL = 3;
gapR = .5;
gapT = .5;
gapB = 1.7;
% gapW = 60;
gapH = 1.5;


subW = (figW-gapL-gapR);
subH = 0.8*subW;

figH = subH+gapT+gapB;

% subplot positioning
subx = gapL/figW;
suby = gapB/figH;

% subx_a = subx_b;
% suby_a = suby_b + (gapH+subH)/figH;

h = figure('name','Overview','units','centimeters','Position', [0, 0, figW, figH])
% set(gcf,'units','centimeters')
subplot('position',[subx suby subW/figW subH/figH])
% plot_planet_profile(col,specs,'a)')
plot_permeability(col,specs,'.')
% subplot('position',[subx_b suby_b subW/figW subH/figH])


%print(gcf, 'Figure_permeability.pdf','-dpdf','-r600');
disp('done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = plot_permeability(col,specs,sub_plot)
% load data
load Freitag1999
load Kawamura2006
load Golden2007
load Ghanbarzadeh2017
load BackGroundPorosity.mat

% Analytic parameterizations
phi = linspace(1e-4,.95,1e2);
d = 1e-3;
k_vBW = d.^2*phi.^2/1600;
k_WW = d^2*phi.^3/200;
k_GHP = d^2*phi.^3/54;

loglog(Freitag1999.phi,Freitag1999.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',6), hold on

%Firn data: Calonne et al. (2022), The Cryosphere
Calonne2022 = readtable("Deff_Keff_data.xlsx");
Calonne2022_rho = table2array(Calonne2022(:,8)); Calonne2022_perm = table2array(Calonne2022(:,17));
loglog(1-Calonne2022_rho./917,Calonne2022_perm,'o','markerfacecolor','w','markeredgecolor',col.green,'MarkerSize',6)

loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
plot(phi,k_GHP,'color',col.blue,'linewidth',specs.linewidth)
plot(phi,d^2*phi.^2.6/180,'color',col.tan,'linewidth',specs.linewidth)
loglog(Freitag1999.phi,Freitag1999.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',6), hold on
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
loglog(Kawamura2006.phi,Kawamura2006.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',6)
loglog(Golden2007.phi,Golden2007.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',6)
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
loglog(Ghanbarzadeh2017.phi30,Ghanbarzadeh2017.k30,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
loglog(Ghanbarzadeh2017.phi60,Ghanbarzadeh2017.k60,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
%patch([NaCl.phi_min NaCl.phi_max NaCl.phi_max NaCl.phi_min],[1e-17 1e-17 1e-9 1e-9],.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.4)
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
loglog(Ghanbarzadeh2017.phi30,Ghanbarzadeh2017.k30,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)
loglog(Ghanbarzadeh2017.phi60,Ghanbarzadeh2017.k60,'o','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',6)


% plot(8.1221e-04*[1 1],[1e-17 1e-9],'k--')
% text(1.7e-3,2.5e-12,'NaCl','fontsize',specs.fontsize_text)
xlim([1e-3 1])
ylim([1e-17 1e-9])
xlabel('$$\phi$$ [-]','fontsize',specs.fontsize_axis_label)
ylabel('$$k [m^2]$$','fontsize',specs.fontsize_axis_label)
legend('Sea ice data','Firn data','Textural Eqbm.','$$\textrm{k}_0 \phi^3, \textrm{k}_0=1.85\times10^{-8} \textrm{m}^2$$','$$\textrm{k}_0 \phi^{2.6}, \textrm{k}_0=5.56\times10^{-9} \textrm{m}^2$$','location','southeast')
text(3e-5,3e-9,sub_plot,'fontsize',specs.fontsize_subpanel)
set(gca,'fontsize',specs.fontsize_axis,'ytick',[1e-17 1e-15 1e-13 1e-11 1e-9])
print(gcf, 'Figure_permeability.pdf','-dpdf','-r600');
end
