%Calculating volume at each time
    set(groot,'defaultAxesFontName','Times')
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultLegendInterpreter','latex')
    set(groot, 'DefaultFigureVisible', 'on');

blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;

d = 10e3;   %Thickness of ice shell [in m]   




load("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeDS = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrDS = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthDS= [Grid.p.yc(pene_ind(1))];

i = 10;

%kc1.85e-8m2
while i<=22200
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrDS = [Time_arrDS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

VolumeDS = [VolumeDS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+10);  %Just staying away from ocean
%pene_depthDS = [pene_depthDS;Grid.p.yc(pene_ind(1))];

[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthDS = [pene_depthDS;pene_depthdummy];

if i<100
    i = i+10;
else
    i = i+100;
end

end

%[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
%Yc_wet = Yc(find(phi>1e-4)); min(Yc_wet(Yc_wet>0))

%Stokes
load("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeS = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrS = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthS= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=100000
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrS = [Time_arrS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeS = [VolumeS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+1);  %Just staying away from ocean
[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthS = [pene_depthS;pene_depthdummy];

i = i+100;
end

%Mid
load("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeM = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrM = [tVec(1)]
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthM= [Grid.p.yc(pene_ind(1))];
i = 2;
while i<=13600
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrM = [Time_arrM;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeM = [VolumeM;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+1);  %Just staying away from ocean
%pene_depthM = [pene_depthM;Grid.p.yc(pene_ind(1))];

[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthM = [pene_depthM;pene_depthdummy];

if i<400
    i = i+1;
else
    i = i+100;
end

end

%Plotting
%Volume plot
figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Darcy-Stokes','Stokes');

%Volume plot
figure()
semilogx(Time_arrDS,pene_depthDS/1e3,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arrS,pene_depthS/1e3,'-','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Penetration depth [km]');
legend('Darcy-Stokes','Stokes');


%Combined
hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'k-','Linewidth',3)
hold on
semilogx(Time_arrM,VolumeM/1e9,'k-.','Linewidth',3)
semilogx(Time_arrS,VolumeS/1e9,'k--','Linewidth',3)
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',3,'color',blue)
semilogx(Time_arrM,VolumeM/1e9,'-.','Linewidth',3,'color',blue)
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',3,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,10 - pene_depthDS*d/1e3,'-','Linewidth',3,'color',red)
hold on
semilogx(Time_arrM,10 - pene_depthM*d/1e3,'-.','Linewidth',3,'color',red)
semilogx(Time_arrS,10 - pene_depthS*d/1e3,'--','Linewidth',3,'color',red)
xlim([1e-4,1e5])
xticks([1e-4,1e-2,1,1e2,1e4])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
%lh = legend('$$\textrm{k}_0=1.85\cdot10^{-8}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-10}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-16}\textrm{m}^2$$','location','southwest');
%set( lh, 'Box', 'off' ) ;
saveas(hh,sprintf('../figures/CombinedVolandpene.png')); 
saveas(hh,sprintf('../figures/CombinedVolandpene.pdf')); 


%Combined Stokes
hh=figure()
yyaxis left
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',blue)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrS,10 - pene_depthS*d/1e3,'-','Linewidth',5,'color',red)
hold on
%ylim([0,45])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
saveas(hh,sprintf('../figures/CombinedStokes.png')); 
saveas(hh,sprintf('../figures/CombinedStokes.pdf')); 

%stop

%Combined Darcy-Stokes
hh=figure()
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',blue)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,10 - pene_depthDS*d/1e3,'-','Linewidth',5,'color',red)
hold on
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]');
saveas(hh,sprintf('../figures/CombinedDarcyStokes.png')); 
saveas(hh,sprintf('../figures/CombinedDarcyStokes.pdf')); 


%Combined Vol
hh=figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',blue)
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('$$\textrm{k}_0=5.6\times10^{-11}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-16}\textrm{m}^2$$','location','southwest');
saveas(hh,sprintf('../figures/CombinedVol.png')); 
saveas(hh,sprintf('../figures/CombinedVol.pdf')); 



%%Combined all k0

%Testing purposes only
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;

%kc1.85e-8m2
while i<=22200
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-8_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-8_temperature_vol.pdf'));

drain_V_array = [phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [1.85e-8]; %in m2

%Testing purposes only kc=1.85e-10 m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;

%kc1.85e-10m2
while i<=13600
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-10_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-10_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-10]; %in m2


%kc1.85e-12m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-12_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=23100
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-12_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-12_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-12_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-12]; %in m2


%kc1.85e-14m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-14_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=15000
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-14_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-14_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-14_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-14]; %in m2

%kc1.85e-16m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=100000
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-16_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-16_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-16]; %in m2


hhh= figure()
semilogx(kc_array, drain_V_array,'ro','Linewidth',2,'color',red,'MarkerSize',12);
xlabel('$$\textrm{k}_0 [\textrm{m}^2]$$');
ylabel('Melt delivered to ocean [km$$^3$$]');
xticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.png')); 
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.pdf'));