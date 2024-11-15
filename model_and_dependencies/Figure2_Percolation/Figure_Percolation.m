% file: Figure1_Overview.m
% author: Marc Hesse
% description: Introductory figure for oxidant percolation paper.
% date: 21 May 2020

clear all,close all, clc
% set(groot,'ScreenPixelsPerInch',218); 
%Analytical 
set(groot,'defaultAxesFontName','Times')
set(groot,'defaultAxesFontSize',20)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultLineLineWidth', 2);
set(0, 'DefaultFigureRenderer', 'painters'); 
col = marc_colors();

col.blue = [57 106 177]./255;

%set(groot,'defaultAxesFontName','Times')
%set(groot,'defaultAxesFontSize',20)
%set(groot,'defaulttextinterpreter','latex')
%set(groot,'defaultAxesTickLabelInterpreter','latex')
%set(groot,'defaultLegendInterpreter','latex')
%set(groot, 'DefaultFigureVisible', 'on');

% Figue size
scsz = get(0,'ScreenSize');
scw = scsz(3); sch = scsz(4);

specs.linewidth = 2.0;
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
load BackGroundPorosity.mat
load Fowler2023.mat

% Analytic parameterizations
phi = linspace(1e-4,.95,1e2);
d = 1e-3;
k_vBW = d.^2*phi.^2/1600;
k_WW = d^2*phi.^3/200;
k_GHP = d^2*phi.^3/54; %original Meyer and Hewitt (2017)

loglog(Freitag1999.phi,Freitag1999.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',12.133), hold on

%Ice data: Fowler 2023
loglog(Fowler2023.fine.phi.mean,Fowler2023.fine.k.mean,'^','markerfacecolor','w','markeredgecolor',col.skyblue,'MarkerSize',6)


%Firn data: Calonne et al. (2022), The Cryosphere
Calonne2022 = readtable("Deff_Keff_data.xlsx");
Calonne2022_rho = table2array(Calonne2022(:,8)); Calonne2022_perm = table2array(Calonne2022(:,17));

%firn data
loglog(1-Calonne2022_rho./917,Calonne2022_perm,'s','markerfacecolor','w','markeredgecolor',col.green,'MarkerSize',12.133)

%simulations
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)

%model
plot(phi,k_GHP,'color',col.black,'linewidth',specs.linewidth)
k_GHPby10 = (1/10)*d^2*phi.^3/54;
k_GHPby1e2= (1/100)*d^2*phi.^3/54;
k_GHPby1e4= (1/10000)*d^2*phi.^3/54;
k_GHPby1e6= (1/1e6)*d^2*phi.^3/54;
plot(phi,k_GHPby10,'color',[col.black,0.75],'linewidth',specs.linewidth)
plot(phi,k_GHPby1e2,'color',[col.black,0.5],'linewidth',specs.linewidth)
plot(phi,k_GHPby1e4,'color',[col.black,0.25],'linewidth',specs.linewidth)
plot(phi,k_GHPby1e6,'color',[col.black,0.1],'linewidth',specs.linewidth)
%plot(phi,d^2*phi.^2.6/180,'color',col.tan,'linewidth',specs.linewidth)

%data
loglog(Freitag1999.phi,Freitag1999.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',12.133), hold on
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Kawamura2006.phi,Kawamura2006.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',12.133)
loglog(Golden2007.phi,Golden2007.k,'o','markerfacecolor','w','markeredgecolor',col.red,'MarkerSize',12.133)
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Ghanbarzadeh2017.phi30,Ghanbarzadeh2017.k30,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Ghanbarzadeh2017.phi60,Ghanbarzadeh2017.k60,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
%patch([NaCl.phi_min NaCl.phi_max NaCl.phi_max NaCl.phi_min],[1e-17 1e-17 1e-9 1e-9],.5*[1 1 1],'EdgeColor','none','FaceAlpha',0.4)
loglog(Ghanbarzadeh2017.phi10,Ghanbarzadeh2017.k10,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Ghanbarzadeh2017.phi30,Ghanbarzadeh2017.k30,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Ghanbarzadeh2017.phi60,Ghanbarzadeh2017.k60,'x','markerfacecolor','w','markeredgecolor',col.blue,'MarkerSize',12.133)
loglog(Fowler2023.coarse.phi.mean,Fowler2023.coarse.k.mean,'^','markerfacecolor','w','markeredgecolor',col.skyblue,'MarkerSize',6)


xlim([1e-3 1])
ylim([1e-17 1e-9])
xlabel('$$\phi$$ [-]','fontsize',specs.fontsize_axis_label,'Interpreter','latex')
ylabel('$$\textrm{k}$$ [m$$^2$$]','fontsize',specs.fontsize_axis_label,'Interpreter','latex')
%legend('Sea ice data','Firn data','Textural Eqbm.','$$\textrm{k}_0 \phi^3, \textrm{k}_0=1.85\times10^{-8} \textrm{m}^2$$','$$\textrm{k}_0 \phi^{2.6}, \textrm{k}_0=5.56\times10^{-9} \textrm{m}^2$$','location','southeast','Interpreter','latex')
lgd = legend('Sea ice data','Poly. ice data','Firn data','Textural Eqbm.','$$\textrm{k}_0=1.85\cdot10^{-8} \textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot10^{-9} \textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot10^{-10} \textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot10^{-12} \textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot10^{-14} \textrm{m}^2$$','location','southeast','Interpreter','latex',BackgroundAlpha=.1)
lgd.NumColumns = 2; lgd.Direction = "reverse"; lgd.IconColumnWidth = 10;
text(3e-5,3e-9,sub_plot,'fontsize',specs.fontsize_subpanel)
set(gca,'fontsize',specs.fontsize_axis,'ytick',[1e-17 1e-15 1e-13 1e-11 1e-9])
print(gcf, 'Figure_permeability.pdf','-dpdf','-r600');
end
