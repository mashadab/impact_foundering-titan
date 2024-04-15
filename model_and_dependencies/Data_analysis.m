%Calculating volume at each time
    set(groot,'defaultAxesFontName','Times')
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultLegendInxterpreter','latex')
    set(groot, 'DefaultFigureVisible', 'on');

blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;

d = 20e3;   %Thickness of ice shell [in m]   


load("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-11_Ea_50_output_100C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
VolumeDS = [];
Time_arrDS = [];
pene_depthDS= [];
i = 100;

%kc5.6e-11m2
while i<=8400
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-11_Ea_50_output_%sC.mat",num2str(i))); %loading file
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
VolumeS = [];
Time_arrS = [];
pene_depthS= [];
i = 100;
while i<=17100
%Individual file
load(sprintf("../Output/Shigeru_impact_Wakita-Stokes_eta0_14kc5.6e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file

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
hh=figure()
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'k-','Linewidth',5)
hold on
semilogx(Time_arrS,VolumeS/1e9,'k--','Linewidth',5)
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',blue)
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,60 - pene_depthDS*d/1e3,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arrS,60 - pene_depthS*d/1e3,'--','Linewidth',5,'color',red)
%ylim([0,45])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
legend('$$\textrm{k}_0=5.6\times10^{-11}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-16}\textrm{m}^2$$','location','southwest');
saveas(hh,sprintf('../figures/CombinedVolandpene.png')); 
saveas(hh,sprintf('../figures/CombinedVolandpene.png.pdf')); 


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

Time_arreminus12 = [Time_arreminus12;tVec(i)] %Time array [years]
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
