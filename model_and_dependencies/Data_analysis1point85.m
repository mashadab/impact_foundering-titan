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

d = 20e3;   %Thickness of ice shell [in m]   




load("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeDS = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrDS = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthDS= [Grid.p.yc(pene_ind(1))];

i = 100;

%kc1.85e-8m2
while i<=62300
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrDS = [Time_arrDS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

VolumeDS = [VolumeDS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthDS = [pene_depthDS;Grid.p.yc(pene_ind(1))];

i = i+100;
end

%Stokes
load("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
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
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrS = [Time_arrS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeS = [VolumeS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];

i = i+100;
end


%Mid
load("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeM = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrM = [tVec(1)]
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthM= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=57400
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrM = [Time_arrM;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeM = [VolumeM;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_depthM = [pene_depthM;Grid.p.yc(pene_ind(1))];

i = i+100;
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
lag = 27;
hh = figure;
hh.Units = 'centimeters';
% [left bottom width height]
hh.Position = [1,1,20,15]; % should be 19 wide
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'k-','Linewidth',3)
hold on
semilogx(Time_arrM,VolumeM/1e9,'k-.','Linewidth',3)
semilogx(Time_arrS(1:end-lag),VolumeS(1:end-lag)/1e9,'k--','Linewidth',3)
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',3,'color',blue)
semilogx(Time_arrM,VolumeM/1e9,'-.','Linewidth',3,'color',blue)
hold on
semilogx(Time_arrS(1:end-lag),VolumeS(1:end-lag)/1e9,'--','Linewidth',3,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,60 - pene_depthDS*d/1e3,'-','Linewidth',3,'color',red)
hold on
semilogx(Time_arrM,60 - pene_depthM*d/1e3,'-.','Linewidth',3,'color',red)
semilogx(Time_arrS(1:end-lag),60 - pene_depthS(1:end-lag)*d/1e3,'--','Linewidth',3,'color',red)
xlim([1e-4,1e5])
xticks([1e-4,1e-2,1e0,1e2,1e4])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
lh = legend('$$\textrm{k}_0=1.85\cdot10^{-8}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-10}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-16}\textrm{m}^2$$','location','southwest');
set( lh, 'Box', 'off' ) ;
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
semilogx(Time_arrS,60 - pene_depthS*d/1e3,'-','Linewidth',5,'color',red)
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
semilogx(Time_arrDS,60 - pene_depthDS*d/1e3,'-','Linewidth',5,'color',red)
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

%k0 = 5.6e-12\textrm{m}^2
Volumeeminus12 = [];
Time_arreminus12 = [];
pene_deptheminus12= [];
i = 100;


while i<=17200
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-12_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arreminus12 = [Time_arreminus12;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
Volumeeminus12 = [Volumeeminus12;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_deptheminus12 = [pene_deptheminus12;Grid.p.yc(pene_ind(1))];

i = i+100
end

%k0 = 5.6e-13\textrm{m}^2
Volumeeminus13 = [];
Time_arreminus13 = [];
pene_deptheminus13= [];
i = 100;

load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-13_Ea_50_output_1C.mat")); %loading file



while i<=38400
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-13_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arreminus13 = [Time_arreminus13;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
Volumeeminus13 = [Volumeeminus13;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_deptheminus13 = [pene_deptheminus13;Grid.p.yc(pene_ind(1))];

i = i+100;
end


%k0 = 5.6e-14\textrm{m}^2
Volumeeminus14 = [];
Time_arreminus14 = [];
pene_deptheminus14= [];
i = 100;
while i<=18000
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-14_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arreminus14 = [Time_arreminus14;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
Volumeeminus14 = [Volumeeminus14;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_deptheminus14 = [pene_deptheminus14;Grid.p.yc(pene_ind(1))];

i = i+100;
end


%Combined Vol
hh=figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arreminus12,Volumeeminus12/1e9,'-','Linewidth',5,'color',brown)
semilogx(Time_arreminus13,Volumeeminus13/1e9,'-','Linewidth',5,'color',green)
semilogx(Time_arreminus14,Volumeeminus14/1e9,'-','Linewidth',5,'color',black)
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',blue)
xlim([5 1e4]);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('$$\textrm{k}_0=5.6\times10^{-11}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-12}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-13}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-14}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-16}\textrm{m}^2$$','location','southwest');
saveas(hh,sprintf('../figures/CombinedVol_morelines.png')); 
saveas(hh,sprintf('../figures/CombinedVol_morelines.pdf'));


%Combined Penetration depth
hh=figure()
semilogx(Time_arrDS,(3*d-pene_depthDS*d)/1e3,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arreminus12,(3*d-pene_deptheminus12*d)/1e3,'-','Linewidth',5,'color',brown)
semilogx(Time_arreminus13,(3*d-pene_deptheminus13*d)/1e3,'-','Linewidth',5,'color',green)
semilogx(Time_arreminus14,(3*d-pene_deptheminus14*d)/1e3,'-','Linewidth',5,'color',black)
semilogx(Time_arrS,(3*d-pene_depthS*d)/1e3,'-','Linewidth',5,'color',blue)
xlim([5 1e4]);
xlabel('Time [years]');
ylabel('Penetration Depth [km]');
legend('$$\textrm{k}_0=5.6\times10^{-11}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-12}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-13}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-14}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-16}\textrm{m}^2$$','location','southwest');
set(gca, 'YDir','reverse')
saveas(hh,sprintf('../figures/CombinedPeneDepth_morelines.png')); 
saveas(hh,sprintf('../figures/CombinedPeneDepth_morelines.pdf'));





%Testing purposes only
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
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
while i<=49700
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
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
while i<=57300
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-12_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=124000
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-12_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-14_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=108300
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-14_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
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
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file
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


drain_V_array_Europa = [25.3205, 23.96, 27.062,26.1107, 26.4404];

lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
semilogx(kc_array, drain_V_array,'ro','Linewidth',2,'color',red,'MarkerSize',12);
hold on
semilogx(kc_array, drain_V_array_Europa,'kx','Linewidth',2,'color',blue,'MarkerSize',12);
legend('Selk crater','Mannann`an crater','location','northeast');
xlabel('$$\textrm{k}_0 [\textrm{m}^2]$$');
ylabel('Melt delivered to ocean [km$$^3$$]');
xticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
ylim([0,35])
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.png')); 
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.pdf'));



%%Time scales for Titan
%Infiltration times
phi_array = linspace(0,1,1000); %Porosity array
kc_array  = logspace(-16,-8,1000); %absolute permeability mid [m2]
[Phi_array, Kc_array] = meshgrid(phi_array,kc_array);
rho_f = 1e3;         %density of water [kg/m3]
mu_f = 1e-3;         %dynamic viscosity [Pa.s]
%grav_Europa =1.315  % gravity: [m/s2] 
grav_Titan = 1.352;   % gravity: [m/s2]
d_Titan = 60e3;       % Titan ice shell thickness [m]    
K_array = Kc_array.*Phi_array.^3*(rho_f - 0)*grav_Titan/(mu_f); 
v_array = K_array./Phi_array; t_inf_array = (d_Titan./v_array)/(3.154e7);


lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
contourf(Phi_array, Kc_array,log10(t_inf_array),20) %kc, phi and infiltration time contour plot
% Show the colorbar
c = colorbar;
c.Label.String = 'Infiltration time [years]';
c.Ticks = -2:16;
c.TickLabels = compose('10^{%d}',c.Ticks);
% Set the log colobar
%set(gca,'ColorScale','log')
set(gca,'yscale','log')
xlabel('Melt fraction, $$\phi$$ [-]');
ylabel('Permeability coefficient $$\textrm{k}_0$$ [m$$^2$$]');
saveas(hhh,sprintf('../figures/infiltration_time_Titan.png')); 
saveas(hhh,sprintf('../figures/infiltration_time_Titan.pdf'));



%%%%Stokes settling at cold temperature
Rgas = 8.314; % universal gas constant, J K^-1 mol^-1
T_b = 273.16; % melting temperature, K
E_a = 50e3; % viscosity activation energy
Apar = E_a/Rgas/T_b; % viscosity exponent
mu_b = 1e14; %basal dynamic viscosity [Pa.s]
T_array =  273.16 + linspace(0,-20,1000);  %Ice shell temperature [K]
volume_array = linspace(50,500,1000);      %[in km3]
[T_array,Volume_array] = meshgrid(T_array,volume_array);
mu_array = mu_b.*exp(Apar*(T_b./T_array-1));%Temperature corrected basal viscosity [Pa.s]
rho_melt  = 1000; %density of the melt (water) [kg/m^3]
rho_matrix= 917;  %density of the bath (ice) [kg/m^3] 
yr2s      = 365.25*24*60*60; %year to second conversion [s/year]
Rad_array    = (3.*Volume_array./(4.*pi)).^(1/3); %[in km]
V_settling = (2/9)*(rho_melt - rho_matrix).*(Rad_array*1e3).^2.*grav_Titan./mu_array.*(yr2s/1e3); %[km/year]
t_settling_without_refreezing = (d_Titan/1e3)./ V_settling; %settling time [years]


%V_settling = (2/9)*(rho_melt - rho_matrix).*(Rad_array*1e3).^2.*grav_Titan./mu.*(yr2s/1e3); %[km/year]
t_settling = d_Titan./ V_settling; %settling time [years]

lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
contourf(T_array-273.16,Volume_array,(t_settling_without_refreezing),20) %kc, phi and infiltration time contour plot
hold on
% Show the colorbar
c = colorbar;
c.Label.String = 'Foundering time [years]';
xlabel('Ice shell temperature [$$^\circ$$C]');
ylabel('Melt volume [km$$^3$$]');
plot(-17,192.97,'ro','Linewidth',2,'color',red,'MarkerSize',12)
saveas(hhh,sprintf('../figures/foundering_time_Titan.png')); 
saveas(hhh,sprintf('../figures/foundering_time_Titan.pdf'));




%%Time scales for Europa
%Infiltration times
phi_array = linspace(0,1,1000); %Porosity array
kc_array  = logspace(-16,-8,1000); %absolute permeability mid [m2]
[Phi_array, Kc_array] = meshgrid(phi_array,kc_array);
rho_f = 1e3;         %density of water [kg/m3]
mu_f = 1e-3;         %dynamic viscosity [Pa.s]
grav_Europa =1.315;  % gravity: [m/s2] 
d_Europa = 10e3;       % Titan ice shell thickness [m]    
K_array = Kc_array.*Phi_array.^3*(rho_f - 0)*grav_Europa/(mu_f); 
v_array = K_array./Phi_array; t_inf_array = (d_Europa./v_array)/(3.154e7);


lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
contourf(Phi_array, Kc_array,log10(t_inf_array),20) %kc, phi and infiltration time contour plot
% Show the colorbar
c = colorbar;
c.Label.String = 'Infiltration time [years]';
c.Ticks = -2:16;
c.TickLabels = compose('10^{%d}',c.Ticks);
% Set the log colobar
%set(gca,'ColorScale','log')
set(gca,'yscale','log')
xlabel('Melt fraction, $$\phi$$ [-]');
ylabel('Permeability coefficient $$\textrm{k}_0$$ [m$$^2$$]');
saveas(hhh,sprintf('../figures/infiltration_time_Europa.png')); 
saveas(hhh,sprintf('../figures/infiltration_time_Europa.pdf'));



%%%%Stokes settling at cold temperature
Rgas = 8.314; % universal gas constant, J K^-1 mol^-1
T_b = 273.16; % melting temperature, K
E_a = 50e3; % viscosity activation energy
Apar = E_a/Rgas/T_b; % viscosity exponent
mu_b = 1e14; %basal dynamic viscosity [Pa.s]
T_array =  273.16 + linspace(0,-20,1000);  %Ice shell temperature [K]
volume_array = linspace(10,100,1000);      %[in km3]
[T_array,Volume_array] = meshgrid(T_array,volume_array);
mu_array = mu_b.*exp(Apar*(T_b./T_array-1));%Temperature corrected basal viscosity [Pa.s]
rho_melt  = 1000; %density of the melt (water) [kg/m^3]
rho_matrix= 917;  %density of the bath (ice) [kg/m^3] 
yr2s      = 365.25*24*60*60; %year to second conversion [s/year]
Rad_array    = (3.*Volume_array./(4.*pi)).^(1/3); %[in km]
V_settling = (2/9)*(rho_melt - rho_matrix).*(Rad_array*1e3).^2.*grav_Europa./mu_array.*(yr2s/1e3); %[km/year]
t_settling_without_refreezing = (d_Europa/1e3)./ V_settling; %settling time [years]


lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
contourf(T_array-273.16,Volume_array,(t_settling_without_refreezing),20) %kc, phi and infiltration time contour plot
hold on
% Show the colorbar
c = colorbar;
c.Label.String = 'Foundering time [years]';
xlabel('Ice shell temperature [$$^\circ$$C]');
ylabel('Melt volume [km$$^3$$]');
plot(0,31,'ro','Linewidth',2,'color',red,'MarkerSize',12)
saveas(hhh,sprintf('../figures/foundering_time_Europa.png')); 
saveas(hhh,sprintf('../figures/foundering_time_Europa.pdf'));