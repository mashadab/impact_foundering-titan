% Change working directory to the folder containing this script
%scriptFullPath = mfilename('fullpath');   % full path of the current script
%[scriptDir, ~, ~] = fileparts(scriptFullPath);  
%cd(scriptDir);


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

%k0 = 1.85e-9\textrm{m}^2
Volumeeminus9 = [];
Time_arreminus9 = [];
pene_deptheminus9= [];
i = 100;


while i<=40000
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-09_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arreminus9 = [Time_arreminus9;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
Volumeeminus9 = [Volumeeminus9;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_deptheminus9 = [pene_deptheminus9;Grid.p.yc(pene_ind(1))];

i = i+100
end


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
while i<=158000%49700
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



%Testing purposes only kc=1.85e-09 m2
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-09_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;

%kc1.85e-09m2
while i<=40000
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-09_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
saveas(hh,sprintf('../figures/1.85e-09_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-09_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-09]; %in m2



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

%kc1.85e-11m2
d = 20e3
load("../Output/Simple-test_eta0_14kc1.85e-11_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=143000
%Individual file
load(sprintf("../Output/Simple-test_eta0_14kc1.85e-11_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+1000;
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
saveas(hh,sprintf('../figures/1.85e-11_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-11_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-11]; %in m2

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


drain_V_array_Titan = [25.3205, 23.96, 27.062,26.1107, 26.4404];

lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide
semilogx(kc_array, drain_V_array,'ro','Linewidth',2,'color',red,'MarkerSize',12);
hold on
semilogx(kc_array, drain_V_array_Titan,'kx','Linewidth',2,'color',blue,'MarkerSize',12);
legend('Selk crater','Mannann`an crater','location','northeast');
xlabel('$$\textrm{k}_0 [\textrm{m}^2]$$');
ylabel('Melt delivered to ocean [km$$^3$$]');
xticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
ylim([0,35])
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.png')); 
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.pdf'));

lag = 27;
hhh = figure;
hhh.Units = 'centimeters';
% [left bottom width height]
hhh.Position = [1,1,20,15]; % should be 19 wide


drain_V_array_Titan_new=         [26.4404, 26.1107, 27.062, 22.679,   23.96,  24.79,   25.3205,   ];
drain_V_array_Titan_uncertainty = [1.3,      1.2,   0.2312, 0.25,    0.2888,   0.2,     0.015 ];
kc_array_Titan = [1.85e-16, 1.85e-14, 1.85e-12, 1.85e-11, 1.85e-10, 1.85e-9, 1.85e-8];

kc_array_Titan = [1.85e-16, 1.85e-14, 1.85e-12, 1.85e-11, 1.85e-10, 1.85e-9, 1.85e-8];
drain_V_array_Titan = [0.8423,0.8077,0.146,0,17.159,32.7964,38.275];
drain_V_array_uncertainty = [2.166,1.8106,2.898,2.33,4.06,1.2,1.23];
%semilogx(kc_array_Titan, drain_V_array_Titan,'ro','Linewidth',2,'color',red,'MarkerSize',12);

ax = axes();
errorbar(ax, kc_array_Titan,drain_V_array_Titan,drain_V_array_uncertainty.*0,drain_V_array_uncertainty,'ro','Linewidth',2,'color',black,'MarkerSize',12);
hold on
errorbar(ax, kc_array_Titan,drain_V_array_Titan_new,drain_V_array_Titan_uncertainty.*0,drain_V_array_Titan_uncertainty,'r^','Linewidth',2,'color',black,'MarkerSize',12);
set(ax, 'XScale', 'log');

legend('Selk crater','Mannann`an crater','location','northwest');
xlabel('Permeability coefficient $$\textrm{k}_0 [\textrm{m}^2]$$');
ylabel('Melt delivered to ocean [km$$^3$$]');
xticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
ylim([0,40])
saveas(hhh,sprintf('../figures/melt_deliveredvsKc_new.png')); 
saveas(hhh,sprintf('../figures/melt_deliveredvsKc_new.pdf'));



%%Time scales for Titan
%Infiltration times
phi_array = linspace(0,1,1000); %Porosity array
kc_array  = logspace(-16,-8,1000); %absolute permeability mid [m2]
[Phi_array, Kc_array] = meshgrid(phi_array,kc_array);
rho_f = 1e3;         %density of water [kg/m3]
mu_f = 1e-3;         %dynamic viscosity [Pa.s]
%grav_Titan =1.315  % gravity: [m/s2] 
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




%%Time scales for Titan
%Infiltration times
phi_array = linspace(0,1,1000); %Porosity array
kc_array  = logspace(-16,-8,1000); %absolute permeability mid [m2]
[Phi_array, Kc_array] = meshgrid(phi_array,kc_array);
rho_f = 1e3;         %density of water [kg/m3]
mu_f = 1e-3;         %dynamic viscosity [Pa.s]
grav_Titan =1.315;  % gravity: [m/s2] 
d_Titan = 10e3;       % Titan ice shell thickness [m]    
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
volume_array = linspace(10,100,1000);      %[in km3]
[T_array,Volume_array] = meshgrid(T_array,volume_array);
mu_array = mu_b.*exp(Apar*(T_b./T_array-1));%Temperature corrected basal viscosity [Pa.s]
rho_melt  = 1000; %density of the melt (water) [kg/m^3]
rho_matrix= 917;  %density of the bath (ice) [kg/m^3] 
yr2s      = 365.25*24*60*60; %year to second conversion [s/year]
Rad_array    = (3.*Volume_array./(4.*pi)).^(1/3); %[in km]
V_settling = (2/9)*(rho_melt - rho_matrix).*(Rad_array*1e3).^2.*grav_Titan./mu_array.*(yr2s/1e3); %[km/year]
t_settling_without_refreezing = (d_Titan/1e3)./ V_settling; %settling time [years]


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
saveas(hhh,sprintf('../figures/foundering_time_Titan.png')); 
saveas(hhh,sprintf('../figures/foundering_time_Titan.pdf'));


%%
%New part
%Combined
tstamp_DS = [1,5000,10000,20000,22500,24000,25000,30000,35000,42000]; %50000
timestamp_DS = []; Volume_timestamp_DS = []; pene_depth_timestamp_DS = [];
tstamp_M = [1,1000,2500,3000,3500,4000,4500,5500,10000,50000];
timestamp_M = []; Volume_timestamp_M = []; pene_depth_timestamp_M = [];
tstamp_S = [1,1500,3000,6000,11000,20000,40000,70000,85000,100000];
timestamp_S = []; Volume_timestamp_S = []; pene_depth_timestamp_S = [];

for j = 1:size(tstamp_DS,2)
    load(sprintf("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(tstamp_DS(j))))
    timestamp_DS = [timestamp_DS,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    if j>1
        if Vol > Volume_timestamp_DS(end)
            Vol = VolumeDS(find(Time_arrDS==tVec(end)));
        end
    end
    Volume_timestamp_DS = [Volume_timestamp_DS,Vol]; %Calculating volume [in m^3]
    
    %Penetration depth
    pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
    pene_ind(pene_ind>interface+100);  %Just staying away from ocean
    pene_depth_timestamp_DS = [pene_depth_timestamp_DS; Grid.p.yc(pene_ind(1))];

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrDS,VolumeDS/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_DS,Volume_timestamp_DS/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)

%semilogx(Time_arrM,VolumeM/1e9,'-.','Linewidth',3,'color',blue)
hold on
%semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',3,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
plot(Time_arrDS,60 - pene_depthDS*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_DS,60 - pene_depth_timestamp_DS*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
%xlim([1e-4,1e5])
%xticks([1e-4,1e-2,1,1e2,1e4])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
%lh = legend('$$\textrm{k}_0=1.85\cdot10^{-8}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-10}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-16}\textrm{m}^2$$','location','southwest');
%set( lh, 'Box', 'off' ) ;
xlim([-0.1,7.5])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newDS.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newDS.pdf')); 

set(gca, 'XScale', 'log');
xlim([1e-3,2e1])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newDS_log.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newDS_log.pdf')); 

%%
%Middle Case
for j = 1:size(tstamp_M,2)
    load(sprintf("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(tstamp_M(j))))
    timestamp_M = [timestamp_M,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    %Penetration depth
    pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
    pene_ind(pene_ind>interface+100);  %Just staying away from ocean

    if j>1
        if Vol > Volume_timestamp_M(end)
            
          Vol = VolumeM(find(Time_arrM==tVec(end)));
          pene_depth_timestamp_M = pene_depthM(find(Time_arrM==tVec(end)));
        end
    end

    Volume_timestamp_M = [Volume_timestamp_M,Vol]; %Calculating volume [in m^3]
    pene_depth_timestamp_M = [pene_depth_timestamp_M; Grid.p.yc(pene_ind(1))];   

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrM,VolumeM/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_M,Volume_timestamp_M/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
plot(Time_arrM,60 - pene_depthM*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_M,60 - pene_depth_timestamp_M*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
xlim([-1,1200])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newM.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newM.pdf')); 

set(gca, 'XScale', 'log');
xlim([1e1,1e4])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newM_log.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newM_log.pdf')); 



%%
%Stokes Case
for j = 1:size(tstamp_S,2)
    load(sprintf("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(tstamp_S(j))))
    timestamp_S = [timestamp_S,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    if j>1
        if Vol > Volume_timestamp_S(end)
            Vol = VolumeS(find(Time_arrS==tVec(end)));
        end
    end
    Volume_timestamp_S = [Volume_timestamp_S,Vol]; %Calculating volume [in m^3]
    
    %Penetration depth
    [Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
    Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
    %pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
    pene_depth_timestamp_S = [pene_depth_timestamp_S; pene_depthdummy];

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrS,VolumeS/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_S,Volume_timestamp_S/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right

%Interp triangle
triang_x = [timestamp_S(6)-50,timestamp_S(4)-50];  % chose location
triang_y = [60 - pene_depth_timestamp_S(6)*d/1e3,60 - pene_depth_timestamp_S(4)*d/1e3];

%aaa1 = plot(triang_x([1,2,2]), triang_y([1,1,2]), 'k-','Linewidth',1.5);
%aaa2 = plot(triang_x, triang_y, 'k-','Linewidth',1.5);

plot(Time_arrS,60 - pene_depthS*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_S,60 - pene_depth_timestamp_S*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
%ylim([3,10.2])
xlim([-100,10000])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newS.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newS.pdf')); 

%delete(aaa1)
%delete(aaa2)
set(gca, 'XScale', 'log');
xlim([1e-1,1e4])
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newS_log.png')); 
saveas(hh,sprintf('../figures/TitanCombinedVolandpene_newS_log.pdf')); 



